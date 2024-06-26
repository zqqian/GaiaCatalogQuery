//
// Created by zqqia on 2021/12/14.
//
//#include "healpix/healpix.cpp"
#include <fstream>
#include <cmath>
#include "algorithm"
#include "CatalogCalMaster.h"

double haversines(double ra1, double dec1,
                  double ra2, double dec2) {
    // cout<<"dec"<<lat1<<"ra"<<lon1<<"dec"<<lat2<<"ra"<<lon2<<endl;
    // distance between latitudes
    // and longitudes

    if (ra1 > 360)ra1 = ra1 - 360;
    if (ra1 < 0)ra1 = ra1 + 360;
    if (ra2 > 360)ra2 = ra2 - 360;
    if (ra2 < 0)ra2 = ra2 + 360;

    double lat1 = dec1;
    double lon1 = ra1;
    double lat2 = dec2;
    double lon2 = ra2;
    double dLat = (lat2 - lat1) *
                  3.14159265358979323846 / 180.0;
    double dLon = (lon2 - lon1) *
                  3.14159265358979323846 / 180.0;

    // convert to radians
    lat1 = (lat1) * 3.14159265358979323846 / 180.0;
    lat2 = (lat2) * 3.14159265358979323846 / 180.0;

    // apply formulae
    double a = pow(sin(dLat / 2), 2) +
               pow(sin(dLon / 2), 2) *
               cos(lat1) * cos(lat2);
    double r = 2 * asin(sqrt(a));//between 0 and 2*pi
    //convert to degree
    double rad = std::fmod(r, 2.0 * 3.14159265358979323846);
    if (rad < 0) {
        rad *= -1;
    }
    double deg = rad * (180. / 3.14159265358979323846);

    return deg;
}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::setSearchRadius(double r) {
    if (r > 0 && r < 360)
        SearchThreshold = r;
}

double ZQQ_AT_CHINA_VO::CatalogCalMaster::getSearchRadius() {
    return SearchThreshold;
}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::build(std::vector<std::shared_ptr<Catalog::line>> BuildTreeVector) {
    min_dis = 99999.0;
    if (!nodelist.empty()) { nodelist.clear(); }

    {
        int id = 0;
        for (auto &lines: BuildTreeVector) {

            wait_to_build_kd k;
            k.coordinate[0] = lines->getRa();
            k.coordinate[1] = lines->getDec();
            k.id = lines->getId();
            k.Line = lines;
            nodelist.push_back(k);
            if (lines->getRa() > 350) {
                k.coordinate[0] = lines->getRa() - 360;
                nodelist.push_back(k);
            }
            if (lines->getRa() < 10) {
                k.coordinate[0] = lines->getRa() + 360;
                nodelist.push_back(k);
            }
            //cout<<k.Line->content<<endl;
        }
    }

    if (!tree_node_list.empty()) { tree_node_list.clear(); }

    tree_node now;
    now.fa = 0;
    tree_node_list.push_back(now);
    tree_build(0, 0, nodelist.size(), 1, 1);
}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::tree_build(int fa, int l, int r, int deep, int type) {
    //cout<<fa<<" "<<l<<" "<<r<<" "<<deep<<endl;
    if (l == r)return;
    int mid = (l + r) >> 1;
    if (l != r - 1)
        std::nth_element(nodelist.begin() + l, nodelist.begin() + mid, nodelist.begin() + r,
                         [=](wait_to_build_kd a, wait_to_build_kd b) -> bool {
                             return a.coordinate[deep % 2] < b.coordinate[deep % 2];
                         });
    tree_node now;
    now.fa = fa;
    now.deep = deep;
    now.is_vaild = 1;
    now.coordinate[0] = nodelist[mid].coordinate[0];
    now.coordinate[1] = nodelist[mid].coordinate[1];
    now.file_id = nodelist[mid].id;
    now.Line = nodelist[mid].Line;
    tree_node_list.push_back(now);
    int id = tree_node_list.size() - 1;
    if (!type) {
        tree_node_list[fa].lson = id;
    } else {
        tree_node_list[fa].rson = id;

    }
    if (l != r - 1) {

        tree_build(id, l, mid, deep + 1, 0);

        tree_build(id, mid + 1, r, deep + 1, 1);

    } else {
        tree_node_list[id].lson = -1;
        tree_node_list[id].rson = -1;

    }
}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::query() {

}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::findAllinThreshold() {
    queryAllCatalogs(ALL_IN_THRESHOLD);
}

void
ZQQ_AT_CHINA_VO::CatalogCalMaster::singleFileQueryAllinThreshold(std::vector<std::shared_ptr<Catalog::line>> l, matched_file &f) {

    quicktimes = 0;
    //  std::vector<std::priority_queue<wait_to_build_kd>> res;
    for (auto  i: l) {
        std::priority_queue<kd_query_result> pq;
        queryKdThreshold(i->getRa(), i->getDec(), 1, pq);
        if (pq.size() > 0) {
            matched_line m;

            m.original = *i;
            while (!pq.empty()) {
                auto p = pq.top();
                matched mm;
                mm.dis = p.dis;
                mm.l = p.Line;
                m.match_details.push_back(mm);
                pq.pop();
            }
            f.line_details.push_back(m);

        } else {
            //cout<<"q:"<<i.getRa()<<" "<<i.getDec()<<" "<<i.content<<"\n";

        }

    }

    return;
}


void ZQQ_AT_CHINA_VO::CatalogCalMaster::queryKdThreshold(double ra, double dec, int now,
                                                         std::priority_queue<kd_query_result> &kp) {
    quicktimes++;
    if (now < 0) { return; }
    if (tree_node_list[now].is_vaild == 0) { return; }
    double cor[2];
    cor[0] = ra;
    cor[1] = dec;
    int closest_id;
    int deep = tree_node_list[now].deep;
    auto &now_node = tree_node_list[now];
    if (cor[deep % 2] > now_node.coordinate[deep % 2]) {
        if (now_node.rson != -1)
        queryKdThreshold(ra, dec, now_node.rson, kp);
    } else {
        if (now_node.lson != -1)
        queryKdThreshold(ra, dec, now_node.lson, kp);

    }

    double dis;
    double min_dis = abs(dec - now_node.coordinate[1]);
    if (min_dis > SearchThreshold) {
        dis = min_dis;

    } else {
        // cout << "haver" << "\n";
        dis = haversines(ra, dec, now_node.coordinate[0], now_node.coordinate[1]);
        havertimes++;
    }
    if (dis < SearchThreshold) {
        kd_query_result k;
        k.id = now_node.file_id;
        k.dis = dis;
        k.Line = now_node.Line;
        kp.push(k);
    }
    if (abs(now_node.coordinate[deep % 2] - cor[deep % 2]) <
        SearchThreshold) {// if distance bigger than radius, it will have match in opposite area probably
        if (cor[deep % 2] < now_node.coordinate[deep % 2]) {
            queryKdThreshold(ra, dec, now_node.rson, kp);
        } else {
            queryKdThreshold(ra, dec, now_node.lson, kp);
        }
    }

}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::printRes() {
    //auto z=&UserVector[0].catalog[0].content;
    //  cout<<*z<<endl;
    int num = 0;
    for (auto f: res) {
        std::cout << "file " << f.filename << "\n";
        for (auto l: f.line_details) {
            if (num > 100)break;
            if (l.match_details.size() == 0)continue;
            std::cout << l.original.content << "\n";
            for (auto m: l.match_details) {
                std::cout << "---\n" << (m.l->content) << "---\n";
                num++;
            }
            std::cout << "------------------------------------------------------\n";
        }
    }

}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::printTot() {


    std::cout << "totInMemory match " << totMatched() << "\n";

}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::findClosest() {
    queryAllCatalogs(CLOSEST);
}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::queryAllCatalogs(match_type type) {

    if (UserVector.size() < 2) {
        std::cout << "file not enough!" << "/n";
        return;
    }

    for (int i = 0; i < UserVector.size(); ++i) {
        build((UserVector[i].getCatalog()));
        matched_file f;

        f.filename = UserVector[i].name;
        for (int j = i + 1; j < UserVector.size(); j++) {
            switch (type) {
                case ALL_IN_THRESHOLD:
                    singleFileQueryAllinThreshold(UserVector[j].getCatalog(), f);
                    break;
                case CLOSEST:
                    singleFileQueryClosest(UserVector[j].getCatalog());
                    break;
                case VIOLENCE:
                    break;
            }

        }
        res.push_back(f);
    }

}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::singleFileQueryClosest(std::vector<std::shared_ptr<Catalog::line>> l) {
    quicktimes = 0;
    //  priority_queue<kd_query_result> pq;

    for (auto i: l) {
        kd_query_result closest;
        closest.dis = 360.0;
        closest.id = -1;
        queryKdClosest(i->getRa(), i->getDec(), 1, closest);
//        if (pq.size() > 0)
//            //cout<<"q:"<<i.getRa()<<" "<<i.getDec()<<" "<<pq.size()<<"\n";
//            if (!pq.empty()) {
//                res.push_back(pq);
//            }


/*output*/
        //   std::cout << i->getRa() << " " << i->getDec() << "\n";
        //   std::cout << closest.dis << " --- " << closest.Line->getRa() << " " << closest.Line->getDec() << "\n";

        if (closest.id != -1 && closest.dis < SearchThreshold)
            //    pq.push(closest);
            ;
        //   Time_cal.push_back(clock());

    }
    // res.push_back(pq);
    return;
}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::queryKdClosest(double ra, double dec, int now, kd_query_result &closest) {
    //cout  <<"now" <<now<< "\n";
    //quicktimes++;
    if (now < 0) { return; }
    if (tree_node_list[now].is_vaild == 0) { return; }
    double cor[2];
    cor[0] = ra;
    cor[1] = dec;
    int closest_id;
    int deep = tree_node_list[now].deep;
    auto &now_node = tree_node_list[now];
    if (cor[deep % 2] > now_node.coordinate[deep % 2]) {
        if (now_node.rson != -1)
        queryKdClosest(ra, dec, now_node.rson, closest);
    } else {
        if (now_node.lson != -1)
        queryKdClosest(ra, dec, now_node.lson, closest);

    }

    double dis;
    //dis = (ra - now_node.coordinate[0]) * (ra - now_node.coordinate[0]) + (dec - now_node.coordinate[1]) * (dec - now_node.coordinate[1]);
    dis = haversines(ra, dec, now_node.coordinate[0], now_node.coordinate[1]);

    if (dis < closest.dis) {
        closest.id = now_node.file_id;
        //for (int i = 0; i < 2; i++)closest.coordinate[i] = now_node.coordinate[i];
        closest.dis = dis;
        closest.Line = now_node.Line;
    }
    if (abs(now_node.coordinate[deep % 2] - cor[deep % 2]) <
        closest.dis) {// if distance bigger than radius, it will have match in opposite area probably
        if (cor[(deep) % 2] > now_node.coordinate[(deep) % 2]) {
            queryKdClosest(ra, dec, now_node.lson, closest);
        } else {
            queryKdClosest(ra, dec, now_node.rson, closest);
        }
    }

}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::findClosestViolence() {
    queryAllCatalogs(VIOLENCE);
}


long long ZQQ_AT_CHINA_VO::CatalogCalMaster::totMatched() {
    long long tot = 0;
    for (auto i: res) {
        for (auto l: i.line_details){
            // totInMemory += l.match_details.size();
            for(auto m:l.match_details){
                if(m.l->content_details[0]!=l.original.content_details[0]){
                    tot++;
                    Log()<<"ERR THERE";
                }
            }
        }


    }
    return tot;
}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::buildAll() {
    min_dis = 99999.0;
    if (!nodelist.empty()) { nodelist.clear(); }
    for (auto &BuildTreeVector: UserVector) {
        int id = 0;
        for (auto &lines: BuildTreeVector.getCatalog()) {
            wait_to_build_kd k;
            k.coordinate[0] = lines->getRa();
            k.coordinate[1] = lines->getDec();
            k.id = id++;
            k.Line=lines;
            nodelist.push_back(k);
            if (lines->getRa() > 350) {
                k.coordinate[0] = lines->getRa() - 360;
                nodelist.push_back(k);
            }
            if (lines->getRa() < 10) {
                k.coordinate[0] = lines->getRa() + 360;
                nodelist.push_back(k);
            }
            //cout<<k.Line->content<<endl;
        }
    }

    if (!tree_node_list.empty()) { tree_node_list.clear(); }

    tree_node now;
    now.fa = 0;
    tree_node_list.push_back(now);
    tree_build(0, 0, nodelist.size(), 1, 1);
}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::singleFileQueryThreshold() {


    //build(UserVector[0].getCatalog());
    matched_file f;

    f.filename = UserVector[0].name;

    singleFileQueryAllinThreshold(UserVector[0].getCatalog(), f);


    res.push_back(f);
}

ZQQ_AT_CHINA_VO::CatalogCalMaster::CatalogCalMaster(ZQQ_AT_CHINA_VO::Catalog c) {
    UserVector.push_back(c);
}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::addCatalog(ZQQ_AT_CHINA_VO::Catalog c) {
    UserVector.push_back(c);
}

void ZQQ_AT_CHINA_VO::CatalogCalMaster::output() {

}

bool ZQQ_AT_CHINA_VO::CatalogCalMaster::queryClosest(double ra, double dec) {
    kd_query_result closest;
    closest.dis = 360.0;
    closest.id = -1;
    queryKdClosest(ra,dec,1,closest);
    //
//    Log()<<ra<<dec<<"closest id:"<<closest.id<<" dis:"<<closest.dis;
//    Log()<<closest.Line->getRa()<<closest.Line->getDec();
//    Log()<<closest.Line->getPix();
    Log()<<"dis"<<closest.dis;
    if(closest.dis>0.01){
//        Log()<<ra<<dec<<"closest id:"<<closest.id<<" dis:"<<closest.dis;
//        Log()<<closest.Line->getRa()<<closest.Line->getDec();
//        Log()<<closest.Line->getPix();
    }
    //
    if(closest.dis<0.01){
        return true;
    }else{

        return false;
    }
}

ZQQ_AT_CHINA_VO::CatalogCalMaster::kd_query_result ZQQ_AT_CHINA_VO::CatalogCalMaster::QueryClosetWithReturn(double  ra,double dec) {
    kd_query_result closest;
    closest.dis = 360.0;
    closest.id = -1;
    queryKdClosest(ra,dec,1,closest);
    return closest;
}

std::vector<ZQQ_AT_CHINA_VO::CatalogCalMaster::kd_query_result>
ZQQ_AT_CHINA_VO::CatalogCalMaster::QueryThresholdWithReturn(double ra, double dec,int lineLimit) {
    std::priority_queue<kd_query_result> pq;
    queryKdThreshold(ra, dec, 1, pq);
    Log()<<"match size"<<pq.size();
    if(pq.size()>0){
        int tot=0;
        std::vector<ZQQ_AT_CHINA_VO::CatalogCalMaster::kd_query_result> k;
        while(!pq.empty()&&tot++<lineLimit){
            k.push_back(pq.top());
            pq.pop();
        }
        return k;
    }
    return std::vector<ZQQ_AT_CHINA_VO::CatalogCalMaster::kd_query_result>{};
}


