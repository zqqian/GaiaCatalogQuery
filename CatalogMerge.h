//
// Created by zqqia on 1/23/2022.
//

#ifndef CROSSMATCHING_CATALOGMERGE_H
#define CROSSMATCHING_CATALOGMERGE_H
#include "CatalogCalMaster.h"
#include "CsvRead.h"

namespace ZQQ_AT_CHINA_VO{
    class CatalogMerge :public CatalogCalMaster{
    public:
        void output()override;
        void merge();
    private:
        void write_to_file(std::string filepath, std::string filename, std::vector<std::shared_ptr<Catalog::line>> content) ;
    };
}



#endif //CROSSMATCHING_CATALOGMERGE_H
