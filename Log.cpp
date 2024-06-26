//
// Created by zqqia on 3/3/2022.
//

#include <iostream>
#include "Log.h"


void ZQQ_AT_CHINA_VO::Log::add(std::string content, ZQQ_AT_CHINA_VO::Log::warningLevel w) {
  //  printConsole(content, w);
    c.append(content+" ");
    level=w;
}

void ZQQ_AT_CHINA_VO::Log::printConsole(std::string content, ZQQ_AT_CHINA_VO::Log::warningLevel w) {
    time_t now = time(0);
    char *dt = ctime(&now);
    auto levelname = magic_enum::enum_name(w);
    std::cout << "[" << levelname << "] : [" << content << "]-->" << dt;
}

ZQQ_AT_CHINA_VO::Log &ZQQ_AT_CHINA_VO::Log::operator<<(std::string c) {
    add(c);
    return *this;
}

ZQQ_AT_CHINA_VO::Log &ZQQ_AT_CHINA_VO::Log::operator<<(int c) {
    add(std::to_string(c));
    return *this;
}

ZQQ_AT_CHINA_VO::Log &ZQQ_AT_CHINA_VO::Log::operator<<(long c) {
    add(std::to_string(c));
    return *this;
}

ZQQ_AT_CHINA_VO::Log &ZQQ_AT_CHINA_VO::Log::operator<<(unsigned long c) {
    add(std::to_string(c));
    return *this;
}

ZQQ_AT_CHINA_VO::Log &ZQQ_AT_CHINA_VO::Log::operator<<(long long int c) {
    add(std::to_string(c));
    return *this;
}

ZQQ_AT_CHINA_VO::Log &ZQQ_AT_CHINA_VO::Log::operator<<(unsigned long long int c) {
    add(std::to_string(c));
    return *this;
}

ZQQ_AT_CHINA_VO::Log &ZQQ_AT_CHINA_VO::Log::operator<<(double c) {
    add(std::to_string(c));
    return *this;
}

ZQQ_AT_CHINA_VO::Log &ZQQ_AT_CHINA_VO::Log::operator<<(float c) {
    add(std::to_string(c));
    return *this;
}

ZQQ_AT_CHINA_VO::Log &ZQQ_AT_CHINA_VO::Log::operator<<(long double c) {
    add(std::to_string(c));
    return *this;
}

ZQQ_AT_CHINA_VO::Log &ZQQ_AT_CHINA_VO::Log::operator<<(unsigned int c) {
    add(std::to_string(c));
    return *this;
}

ZQQ_AT_CHINA_VO::Log &ZQQ_AT_CHINA_VO::Log::operator<<(ZQQ_AT_CHINA_VO::Log::warningLevel w) {
    level=w;
    return *this;
}


