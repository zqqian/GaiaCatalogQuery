
#include "ConsoleRead.h"
void ZQQ_AT_CHINA_VO::ConsoleRead::input(std::string name) {
    int n;
    std::cin >> n;

    for (int i = 0; i < n; ++i) {
        double ra,dec;
        std::cin>>ra>>dec;
        std::shared_ptr<line> l;
        l->setRa(ra);
        l->setDec(dec);
        l->content="111";
        this->addLine(l);
    }
    Catalog d;
    double ra,dec;
    std::cin>>ra>>dec;
    d.addLine(0,ra,dec,"222");


}

void ZQQ_AT_CHINA_VO::ConsoleRead::output() {
    //printRes();
   // printTot();
}

