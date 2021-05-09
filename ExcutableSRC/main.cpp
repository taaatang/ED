#include "Measure/measure.hpp"

int main(int argc, char* argv[]) {
    int workerID{0}, workerNum{1};
    bool isMaster;
    init(workerID, workerNum, isMaster);
    Timer timer(isMaster);

    assert_msg(argc > 1, "No inputDir!");
    std::string inputDir(argv[1]);
    std::vector<std::string> jobs;
    // cmd line jobs
    for (int i = 2; i < argc; ++i) {
        jobs.push_back(std::string(argv[i]));
    }

    System<cdouble> sys(inputDir, isMaster);
    sys.diag();
    // input file jobs
    if (jobs.empty()) { 
        for (auto job : sys.measurePara.getvecs("all")) {
            if (sys.measure(job)) {
                jobs.push_back(job);
            }
        }
    }

    for (auto job : jobs) {
        if (isMaster) {
            printLine();
            std::cout << "Begin " << job << '\n';
        }
        timer.tik();
        compute(sys, job, workerID, workerNum);
        timer.print(job);
    }

    MPI_Finalize();
    return 0;
}