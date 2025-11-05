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

std::string version = "0.0.7";

//global
std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); // immediately start the clock when the program is run

int verbose_flag;
Log lg;
std::vector<Log> logs;
int tabular_flag;
std::mutex mtx;
ThreadPool<std::function<bool()>> threadPool;

int main(int argc, char **argv) {
    
    short int c; // optarg
    short unsigned int pos_op = 1; // optional arguments
    bool arguments = true;
    bool isPipe = false; // to check if input is from pipe
    UserInputRdeval userInput; // initialize input object
    std::string cmd;
    char helpStr[] = "rdeval input.fa*[.gz]|bam|cram|rd [expected genome size as uint]";
    
    if (argc == 1) { // gfastats with no arguments
        printf("%s\n", helpStr);
        exit(0);
    }
    static struct option long_options[] = { // struct mapping long options
        {"homopolymer-compress", required_argument, 0, 0},
        {"sample", required_argument, 0, 0},
        {"random-seed", required_argument, 0, 0},
		{"parallel-files", required_argument, 0, 0},
        {"decompression-threads", required_argument, 0, 0},
        {"compression-threads", required_argument, 0, 0},
        {"sequence-report",no_argument, 0, 0},
        
        {"exclude-list", required_argument, 0, 'e'},
        {"filter", required_argument, 0, 'f'},
        {"include-list", required_argument, 0, 'i'},
        {"threads", required_argument, 0, 'j'},
		{"max-memory", required_argument, 0, 'm'},
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
        c = getopt_long(argc, argv, "-:e:f:m:i:j:o:r:s:q:c:vh",
                        long_options, &option_index);
        if (c == -1) // exit the loop if run out of options
            break;
        switch (c) {
            case ':': // handle options without arguments
                switch (optopt) { // the command line option last matched
                    default:
                        fprintf(stderr, "option -%c is missing a required argument\n", optopt);
                        return EXIT_FAILURE;
                }
                break;
            case 0: // case for long options without short options
                if(strcmp(long_options[option_index].name,"homopolymer-compress") == 0)
                    userInput.hc_cutoff = atoi(optarg);
                if(strcmp(long_options[option_index].name,"sample") == 0) {
                    try {
                        userInput.ratio = std::stof(optarg);
                    } catch (const std::exception& e) {
                        std::cerr << "Invalid float value: " << optarg << std::endl;
                        return EXIT_FAILURE;
                    }
                    if (userInput.ratio > 1 || userInput.ratio < 0)
                        fprintf(stderr, "option -%c needs a float between 0 and 1\n", optopt);
                }
                if(strcmp(long_options[option_index].name,"random-seed") == 0)
                    userInput.randSeed = atoi(optarg);
                if(strcmp(long_options[option_index].name,"sequence-report") == 0) {
                    userInput.content_flag = 1;
                    userInput.stats_flag = false;
                    userInput.outSize_flag = false;
                    userInput.quality_flag = false;
                }
				if(strcmp(long_options[option_index].name,"parallel-files") == 0)
					userInput.parallel_files = atoi(optarg);
                if(strcmp(long_options[option_index].name,"decompression-threads") == 0)
                    userInput.decompression_threads = atoi(optarg);
                if(strcmp(long_options[option_index].name,"compression-threads") == 0)
                    userInput.compression_threads = atoi(optarg);
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
                }
                break;
            case 'e': // list exclude
                ifFileExists(optarg);
                userInput.inBedExclude = optarg;
                break;
            case 'f' : // filtering input
                userInput.filter = optarg;
                rmChrFromStr(userInput.filter, "'\\");
                break;
            case 'i': // list include
                ifFileExists(optarg);
                userInput.inBedInclude = optarg;
                break;
			case 'm': // max memory
				userInput.maxMem = atoi(optarg);
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
            case 'o':
                userInput.outFiles.push_back(optarg);
                break;
            case 'j': // max threads
                userInput.maxThreads = atoi(optarg);
                break;
            case 'v': // software version
                printf("rdeval v%s\n", version.c_str());
                printf("Giulio Formenti giulio.formenti@gmail.com\n");
                exit(0);
            case 'h': // help
                printf("%s", helpStr);
                printf("\nOptions:\n");
                printf("\t--sequence-report generates a per-read report\n");
                printf("\t-e --exclude-list <file> generates output on a excluding list of headers.\n");
                printf("\t-f --filter <exp> filter reads using <exp> in quotes, e.g. 'l>10' for longer than 10bp or 'l>10 & q>10' to further exclude reads by quality (default: none).\n");
                printf("\t-i --include-list <file> generates output on a subset list of headers.\n");
                printf("\t-o --out-format <file> output file (fa*[.gz], bam, cram, rd). Optionally write reads to file or generate rd summary file.\n");
                printf("\t-q --quality q|a generates list of average quality for each read (q) or both length and quality (a).\n");
                printf("\t-r --input-reads <file1> <file2> <file n> input file (fa*[.gz], bam, cram, rd).\n");
                printf("\t-s --out-size u|s|h|c  generates size list (unsorted|sorted|histogram|inverse cumulative table).\n");
                printf("\t--homopolymer-compress <int> compress all the homopolymers longer than n in the input.\n");
                printf("\t--sample <float> fraction of reads to subsample.\n");
                printf("\t--random-seed <int> an optional random seed to make subsampling reproducible.\n");
                printf("\t--md5 print md5 of .rd files.\n");
				printf("\t--parallel-files <int> numbers of files that can be opened and processed in parallel (default:4).\n");
				printf("\t--decompression-threads <int> numbers of decompression threads used by htslib for bam/cram (default:4).\n");
				printf("\t--compression-threads <int> numbers of compression threads used by htslib for bam/cram (default:6).\n");
                printf("\t--tabular tabular output.\n");
                printf("\t--verbose verbose output.\n");
				printf("\t-m --max-memory <int> max number of bases in ring buffer (default:1000000).\n");
				printf("\t-j --threads <int> numbers of threads (default:8).\n");

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
	maxMem = (userInput.maxMem == 0 ? get_mem_total(3) * 0.9 : userInput.maxMem); // set memory limit
    in.read();
    threadPool.join(); // join threads
    exit(EXIT_SUCCESS);
}
