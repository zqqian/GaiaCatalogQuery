#pragma once
#ifndef CONSOLE_READ_H
#define CONSOLE_READ_H
#include "CatalogCalMaster.h"

namespace ZQQ_AT_CHINA_VO {
    class ConsoleRead : public Catalog {
    public:
        void input(std::string name);
        void output()override;
    };
}

#endif