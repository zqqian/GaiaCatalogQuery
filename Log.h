//
// Created by zqqia on 3/3/2022.
//

#ifndef CROSSMATCHING_LOG_H
#define CROSSMATCHING_LOG_H

#include <string>
#include "magic_enum.h"

namespace ZQQ_AT_CHINA_VO {
    class Log {
    public:
        enum warningLevel {
            DEBUG,
            INFO,
            NOTIFY,
            WARNING,
            ERROR,
            FATAL
        };
        warningLevel screenPrintLevel;
        warningLevel level;

        Log() {

        #ifndef NDEBUG
            screenPrintLevel = DEBUG;
        #else
            screenPrintLevel=DEBUG;
        #endif
        }

        ~Log() {
            if (level >= screenPrintLevel) {
                printConsole(c, level);
            }
            if (level >= FATAL) {
                exit(-1);
            }
            if (level >= ERROR) {
                errno=-1;
                perror("ERROR OCCURRED! \nPRESS ANY CHARACTER TO CONTINUE");
                int c;
                std::cin >> c;
            }


        }


        std::string c;

        void add(std::string content, warningLevel w = NOTIFY);

        Log &operator<<(std::string c);

        Log &operator<<(warningLevel w);

        Log &operator<<(int c);

        Log &operator<<(unsigned int c);

        Log &operator<<(long c);

        Log &operator<<(unsigned long c);

        Log &operator<<(long long c);

        Log &operator<<(unsigned long long c);

        Log &operator<<(double c);

        Log &operator<<(float c);

        Log &operator<<(long double c);


    private:
        void printConsole(std::string content, warningLevel w);
    };
}


#endif //CROSSMATCHING_LOG_H
