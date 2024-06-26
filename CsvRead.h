//
// Created by zqqia on 2021/12/23.
//
#pragma once

#ifndef CROSSMATCHING_CSVREAD_H
#define CROSSMATCHING_CSVREAD_H

#include <string.h>
#include "CatalogCalMaster.h"
#include "algorithm"

namespace ZQQ_AT_CHINA_VO {
/*
 *WARNING! THE CONTENT NUMBER IS START FROM  0 !
 * */
    class CSVRow {
    public:
        bool LineParseMode = {true};
        char seperater = ',';

        CSVRow(bool lineParseMode) : LineParseMode(lineParseMode) {}

        CSVRow() {}

        std::string_view operator[](std::size_t index) const {

            return std::string_view(&m_line[m_data_begin[index]], m_data_end[index] - m_data_begin[index] + 1);

        }

        std::size_t size() const {
            return m_data_begin.size();
        }

        std::string rawString() {
            return m_line;
        }

        void set_string(std::string s) {
            m_line = s;
        }

        virtual bool set_m_data() {
            return true;
        }

        bool readNextRow(FILE *fp) {
            char str[10000];
            auto res = fgets(str, 10000, fp);
            if (res == nullptr)return false;
            int slen = strlen(str);
            //  std::cout<<slen<<std::endl;
            if (str[slen - 1] == '\n') {
                if (str[slen - 2] == '\r') {
                    str[slen - 2] = '\0';
                }
                str[slen - 1] = '\0';
            }

            m_line = str;
            if (LineParseMode)
                return readString(m_line);
            else
                return true;
        }

        void setM_data_begin(const std::vector<int> &s) {
            m_data_begin.assign(s.begin(), s.end());
        }

        void setM_data_end(const std::vector<int> &s) {
            m_data_end.assign(s.begin(), s.end());
        }

        bool readString(std::string &str) {

            // std::getline(str, m_line);
            //zqq
            m_data_begin.clear();
            m_data_end.clear();
            std::string::size_type pos = 0;
            while (m_line[pos] == ' ')pos++;// in case of many spaces in the begin of line
            m_data_begin.emplace_back(pos);
            while ((pos = m_line.find(seperater, pos)) != std::string::npos) {
                m_data_end.emplace_back(pos - 1);
                ++pos;
                while (m_line[pos] == ' ')pos++;
                m_data_begin.emplace_back(pos);
            }
            pos = m_line.size() - 1;
            while (m_line[pos] == ' ')pos--;
            m_data_end.emplace_back(pos);
            return true;
        }

        // the number of data column
        int dataSize() {
            return m_data_begin.size();
        };

    private:
        std::string m_line;
        std::vector<int> m_data_begin;
        std::vector<int> m_data_end;
    };


    class CsvRead : public Catalog {
    public:
        CsvRead();

        CsvRead(std::string Catalogname);

        CsvRead(std::string Catalogname, int headerLine, int lineLimit);

        CsvRead(std::string Catalogname, double epoch);

        CsvRead(std::string Catalogname, std::string strindexType, long long MEMORY_SIZE_LIMIT);

        static const int NO_HEADER = 0;
        static const int FIRST_LINE_HEADER = 1;
        static const int unlimitLineNumber = -1;
        bool SAVE_TO_CATALOG = true;
        enum {
            READ_FROM_HEADER,
            SET_MANUALLY,
            DEFAULT
        } RaDecType = READ_FROM_HEADER;
        bool HaveProperMotion;

        void input(std::string name);

        void inputSingleCsv(std::string name, int headerLine = FIRST_LINE_HEADER, int lineLimit = -1);

        virtual bool lineDeal(Catalog::line &l);

        void inputdir(std::string name, std::string suffix = "csv");

        void output();

        virtual void setRaLocDefault();

        void setRcDecLocFromHeader(std::vector<std::string> header);

        void setRaDecNameInHeader(std::string ra_name, std::string dec_name);

        void readCatalogSingleFileFromTemp(int uniq);

        std::vector<std::string> headerDeal();

        virtual bool RaDecDeal(Catalog::line &l);

        void PmDeal(Catalog::line &l);

        CSVRow row;

        void mergeCatalog(std::string dir);

        int wirteCatalog(int nowlevel, int ipix, std::string Dir) override;

        void write_to_file(std::string filepath, std::string filename, std::string content) {
            // auto startTime = std::chrono::high_resolution_clock::now();
            if (filepath.size() == 0)return;
            if (filepath[filepath.size() - 1] != '/') {
                filepath.append("/");
            }
#if defined(_WIN32) || defined(_WIN64)
            if (_access(filepath.c_str(), 00) == -1) {
                _mkdir(filepath.c_str());
            }
#endif
#if   defined(__linux__)
            if (access(filepath.c_str(), 0) == -1) {
            mkdir(filepath.c_str(),S_IRWXU);
        }

#endif
            std::ofstream outFile;
            const size_t bufsize = 256 * 1024;
            char buf[bufsize];
            outFile.rdbuf()->pubsetbuf(buf, bufsize);
            outFile.open((filepath + filename).c_str(), std::ios_base::app);
            outFile.write(content.c_str(), content.size());

            outFile.close();
            //auto endTime = std::chrono::high_resolution_clock::now();

            Log() << "Write to " << (filepath + filename) << " size " << content.size()
                  << "\n";// " cost time "<< std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " ms \n";
        }

        static std::vector<std::string> get_file_list(std::string name, std::string suffix = "csv");

        std::string temp_dir;
        std::string final_split_dir;
        int64_t CatalogTotRead = 0;

        int GaiaPbfTimes = 0;
    protected:
        int raLoc{};
        int decLoc{};
        struct {
            double pmrLoc = -1, pmdLoc = -1, pxLoc = -1, rvLoc = -1;
        } pmLoc;

        bool checkLocLegal(int loc);

        std::vector<int> GAIApixList;
    private:
        virtual void writeGAIA();

        std::string raname;
        std::string decname;

    };


}


#endif //CROSSMATCHING_CSVREAD_H
