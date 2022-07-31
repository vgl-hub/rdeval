#include <main.h>

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
    
    UserInput userInput; // initialize input object
    
    std::string cmd;
    
    if (argc == 1) { // gfastats with no arguments
            
        printf("mytool [command]\n-h for additional help.\n");
        exit(0);
        
    }
    
    static struct option long_options[] = { // struct mapping long options
        {"input-reads", required_argument, 0, 'r'},
        
        {"cmd", no_argument, &cmd_flag, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        
        {0, 0, 0, 0}
    };
    
    while (arguments) { // loop through argv
        
        int option_index = 0;
        
        c = getopt_long(argc, argv, "-:r:vh",
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
            default: // handle positional arguments
                
            case 0: // case for long options without short options
                
//                if (strcmp(long_options[option_index].name,"line-length") == 0)
//                  splitLength = atoi(optarg);
                
                break;
                
            case 'r': // input reads
                
                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    userInput.pipeType = 'r'; // pipe input is a sequence
                
                }else{ // input is a regular file
                    
                    ifFileExists(optarg);
                    userInput.iReadFileArg = optarg;
                    stats_flag = true;
                    
                }
                    
                break;
                
            case 'v': // software version
                printf("mytool v%s\n", version.c_str());
                printf("Giulio Formenti giulio.formenti@gmail.com\n");
                exit(0);
                
            case 'h': // help
                printf("mytool [command]\n");
                printf("\nOptions:\n");
                printf("-r --reads <file> input file (fasta, fastq [.gz]). Optional reads. Summary statistics will be generated.\n");
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
        
        inReads.report();
        
    }
    
    threadPool.join(); // join threads
    
    exit(EXIT_SUCCESS);
    
}
