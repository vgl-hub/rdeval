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

std::string version = "0.0.2";

//global
std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); // immediately start the clock when the program is run

int verbose_flag;
Log lg;
std::vector<Log> logs;
int tabular_flag;

int maxThreads = 3;
std::mutex mtx;
ThreadPool<std::function<bool()>> threadPool;

int main(int argc, char **argv) {
    
    short int c; // optarg
    short unsigned int pos_op = 1; // optional arguments
    
    bool arguments = true;
    bool isPipe = false; // to check if input is from pipe
    
    UserInputRdeval userInput; // initialize input object
    
    std::string cmd;

    char helpStr[] = "rdeval input.fa*[.gz]|bam|cram [expected genome size]";
    
    if (argc == 1) { // gfastats with no arguments
            
        printf("%s\n", helpStr);
        exit(0);
        
    }
    
    static struct option long_options[] = { // struct mapping long options
        {"homopolymer-compress", required_argument, 0, 0},
        {"content",required_argument, 0, 'c'},
        {"filter", required_argument, 0, 'f'},
        {"threads", required_argument, 0, 'j'},
        {"out-format", required_argument, 0, 'o'},
        {"quality", required_argument, 0, 'q'},
        {"input-reads", required_argument, 0, 'r'},
        {"out-size", required_argument, 0, 's'},
        {"md5", no_argument, &userInput.md5_flag, 1},
        {"tabular", no_argument, &tabular_flag, 1},
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
                
                if(strcmp(long_options[option_index].name,"homopolymer-compress") == 0)
                    userInput.hc_cutoff = atoi(optarg);
                
                break;
            
            default: // handle positional arguments
                
                if (isInt(optarg)) { // if the positional argument is a number, it is likely the expected genome size
                    userInput.gSize = atoll(optarg); pos_op++;
                    break;
                }
                /* fall through */
                
            case 'r': // input files
                
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
                printf("\t-j --threads <n> numbers of threads (default:3).\n");
                printf("\t-f --filter <exp> filter reads using <exp> in quotes, e.g. 'l>10' for longer than 10bp or 'l>10 & q>10' to further exclude reads by quality (default: none).\n");
                printf("\t-c --content a|g|t|n generates a list of sequences and their ATCGN base content; all bases, GC content, AT content, Ns (default:a).\n");
                printf("\t-o --out-format <file> output file (fa*[.gz], rd). Optionally write reads to file or generate rd summary file.\n");
                printf("\t-q --quality c|l a generates list of average quality for each read (c) or both length and quality (c).\n");
                printf("\t-r --reads <file1> <file2> <file n> input file (fasta|fastq [.gz], bam, cram).\n");
                printf("\t-s --out-size u|s|h|c  generates size list (unsorted|sorted|histogram|inverse cumulative table).\n");
                printf("\t--homopolymer-compress <n> compress all the homopolymers longer than n in the input.\n");
                printf("\t--md5 print md5 of .rd files.\n");
                printf("\t--tabular tabular output.\n");
                printf("\t--verbose verbose output.\n");
                printf("\t-v --version software version.\n");
                printf("\t--cmd print $0 to stdout.\n");
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
