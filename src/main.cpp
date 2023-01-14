#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <vector>
#include <getopt.h>
#include <fstream>

#include "global.h"
#include "log.h"
#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "gfa-lines.h"
#include "reads.h"
#include "stream-obj.h"
#include "input.h"

std::string version = "0.0.1";

//global
std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); // immediately start the clock when the program is run

int hc_flag;
int hc_cutoff;
short int tabular_flag;
int cmd_flag;
int verbose_flag;
int outBubbles_flag;
int stats_flag;
int outSize_flag;
int filterInput = 0;
int quality_flag;
int maxThreads = 0;
int discoverPaths_flag;

std::mutex mtx;
ThreadPool<std::function<bool()>> threadPool;
Log lg;

int main(int argc, char **argv) {
    
    short int c; // optarg
    short unsigned int pos_op = 1; // optional arguments
    
    bool arguments = true;
    bool isPipe = false; // to check if input is from pipe
    
    unsigned long long int gSize = 0; // expected genome size, with 0 NG/LG* statistics are not computed
    
    UserInputRdeval userInput; // initialize input object

    char sizeOutType = 'u'; //default output from this flag is unsorted sizes 
    char qualityOut = 'a'; // average quality per read 
    
    std::string cmd;

    char helpStr[] = "rdeval input.[fasta|fastq|gfa][.gz] [expected genome size]";
    
    if (argc == 1) { // gfastats with no arguments
            
        printf("%s", helpStr);
        exit(0);
        
    }
    
    static struct option long_options[] = { // struct mapping long options
        {"input-reads", required_argument, 0, 'r'},
        
        {"threads", required_argument, 0, 'j'},

        {"filter", required_argument, 0, 'f'},
        {"out-size", required_argument, 0, 's'},

        {"quality", required_argument, 0, 'q'},
        
        {"verbose", no_argument, &verbose_flag, 1},
        {"cmd", no_argument, &cmd_flag, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        
        {0, 0, 0, 0}
    };
    
    while (arguments) { // loop through argv
        
        int option_index = 0;
        
        c = getopt_long(argc, argv, "-:r:j:f:s:q:vh",
                        long_options, &option_index);
        
        if (c == -1) { // exit the loop if run out of options
            break;
            
        }

        switch (c) {
            case ':': // handle options without arguments
                switch (optopt) { // the command line option last matched
                    case 'b':
                        break;
                        
                    default:
                        fprintf(stderr, "option -%c is missing a required argument\n", optopt);
                        return EXIT_FAILURE;
                }
                break;
                
            case 0: // case for long options without short options
                
//                if (strcmp(long_options[option_index].name,"line-length") == 0)
//                  splitLength = atoi(optarg);
                
                break;
            
            default: // handle positional arguments
                
                if (isInt(optarg)) { // if the positional argument is a number, it is likely the expected genome size
                    
                    gSize = atoll(optarg); pos_op++;
                    
                    break; 
                    
                }
                /* fall through */
                
            case 'r': // input reads
                
                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    userInput.pipeType = 'r'; // pipe input is a sequence
                
                }else{ // input is a regular file
                    
                    optind--;
                    for( ;optind < argc && *argv[optind] != '-' && !isInt(argv[optind]); optind++){
                        
                        ifFileExists(argv[optind]);
                        userInput.iReadFileArg.push_back(argv[optind]);
                        
                    }
                    
                    stats_flag = true;
                    
                }
                    
                break;
                
            case 'j': // max threads
                maxThreads = atoi(optarg);
                stats_flag = 1;
                break;

            case 'f' : //filtering input
                
                userInput.filter = optarg;
                
                userInput.filter.erase(remove(userInput.filter.begin(), userInput.filter.end(), '\''));
                userInput.filter.erase(remove(userInput.filter.begin(), userInput.filter.end(), '\\'));

                if (!((userInput.filter[0] == '>' || userInput.filter[0] == '<' || userInput.filter[0] == '=') && isInt(userInput.filter.substr(1)))) {
                    printf ("Could not parse filter: %s \n", userInput.filter.c_str());
                    exit(0);
                }
                break;

            case 's':
                sizeOutType = *optarg;
                outSize_flag = 1;
                stats_flag = false;
                quality_flag = false;
                break;

            case 'q':
                qualityOut = *optarg;
                quality_flag = 1;
                stats_flag = false;
                outSize_flag = false;
                break;


                
            case 'v': // software version
                printf("rdeval v%s\n", version.c_str());
                printf("Giulio Formenti giulio.formenti@gmail.com\n");
                exit(0);
                
            case 'h': // help
                printf("%s", helpStr);
                printf("\nOptions:\n");
                printf("-j --threads <n> numbers of threads (default:max).\n");
                printf("-f --filter <n> minimum length for retention (default:0).\n");
                printf("-s --out-size u|s|h|c  generates size list (unsorted|sorted|histogram|inverse cummulative table).\n");
                printf("-q --quality a generates list of average quality for each read.\n");
                printf("-r --reads <file1> <file2> <file n> input file (fasta, fastq [.gz]). Optional reads. Summary statistics will be generated.\n");
                printf("--verbose verbose output.\n");
                printf("-v --version software version.\n");
                printf("--cmd print $0 to stdout.\n");
                exit(0);
        }
        
        if    (argc == 2 || // handle various cases in which the output should include summary stats
              (argc == 3 && pos_op == 2) ||
              (argc == 4 && pos_op == 3)) {
            
        }
        
    }
    
    if (cmd_flag) { // print command line
        for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
        
    }
    Input in;
    in.load(userInput); // load user input

    lg.verbose("Loaded user input");
    
    InReads inReads; // initialize sequence collection object
    lg.verbose("Read object generated");
    
    threadPool.init(maxThreads); // initialize threadpool

    in.read(inReads); // read input content to inReads container
    
    jobWait(threadPool);

    if (stats_flag) { // output summary statistics
        

        inReads.report(gSize);
        
    }

    else if (outSize_flag) {
        
        inReads.printReadLengths(sizeOutType);
    }
    else if (quality_flag) {

        inReads.printQualities(qualityOut);

    }
    
    threadPool.join(); // join threads
    
    exit(EXIT_SUCCESS);
    
}
