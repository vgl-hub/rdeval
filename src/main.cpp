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

int verbose_flag;
Log lg;
std::vector<Log> logs;
int tabular_flag;

int maxThreads = 0;
std::mutex mtx;
ThreadPool<std::function<bool()>> threadPool;

int main(int argc, char **argv) {
    
    short int c; // optarg
    short unsigned int pos_op = 1; // optional arguments
    
    bool arguments = true;
    bool isPipe = false; // to check if input is from pipe
    
    UserInputRdeval userInput; // initialize input object
    
    std::string cmd;

    char helpStr[] = "rdeval input.[fasta|fastq|gfa][.gz] [expected genome size]\n";
    
    if (argc == 1) { // gfastats with no arguments
            
        printf("%s", helpStr);
        exit(0);
        
    }
    
    static struct option long_options[] = { // struct mapping long options
        {"content",required_argument, 0, 'c'},
        {"filter", required_argument, 0, 'f'},
        {"threads", required_argument, 0, 'j'},
        {"out-format", required_argument, 0, 'o'},
        {"quality", required_argument, 0, 'q'},
        {"input-reads", required_argument, 0, 'r'},
        {"out-size", required_argument, 0, 's'},
        {"verbose", no_argument, &verbose_flag, 1},
        {"cmd", no_argument, &userInput.cmd_flag, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        
        {0, 0, 0, 0}
    };
    
    while (arguments) { // loop through argv
        
        int option_index = 0;
        
        c = getopt_long(argc, argv, "-:r:j:f:o:s:q:c:vh",
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
                    userInput.gSize = atoll(optarg); pos_op++;
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
                        userInput.inFiles.push_back(argv[optind]);
                        
                    }
                    
                    userInput.stats_flag = true;
                    
                }
                    
                break;
                
            case 'j': // max threads
                maxThreads = atoi(optarg);
                userInput.stats_flag = 1;
                break;

            case 'f' : //filtering input
                
                userInput.filter = optarg;
                
                rmChrFromStr(userInput.filter, "'\\");

                if (!(userInput.filter[0] == '>' || userInput.filter[0] == '<' || userInput.filter[0] == '=') || !isInt(userInput.filter.substr(1))) {
                    printf ("Could not parse filter: %s \n", userInput.filter.c_str());
                    exit(0);
                }
                break;

            case 's':
                userInput.sizeOutType = *optarg;
                userInput.outSize_flag = 1;
                userInput.stats_flag = false;
                userInput.quality_flag = false;
                break;

            case 'q':
                userInput.qualityOut = *optarg;
                userInput.quality_flag = 1;
                userInput.stats_flag = false;
                userInput.outSize_flag = false;
                break;

            case 'c':
                userInput.content = *optarg;
                userInput.content_flag = 1;
                userInput.stats_flag = false;
                userInput.outSize_flag = false;
                userInput.quality_flag = false;
                break;

            case 'o':
                userInput.outFiles.push_back(optarg);
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
                printf("-c --content a|g|t|n generates a list of sequences and their ATCGN base content; all bases, GC content, AT content, Ns (default:a).\n");
                printf("-o --out-format <file> output file (fasta, fastq [.gz]). Optionally write reads to file.\n");
                printf("-q --quality a generates list of average quality for each read.\n");
                printf("-r --reads <file1> <file2> <file n> input file (fasta, fastq [.gz]).\n");
                printf("-s --out-size u|s|h|c  generates size list (unsorted|sorted|histogram|inverse cummulative table).\n");
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
    
    if (userInput.cmd_flag) { // print command line
        for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
    }
    Input in;
    in.load(userInput); // load user input
    lg.verbose("Loaded user input");
    in.read();
    threadPool.join(); // join threads
    exit(EXIT_SUCCESS);
}
