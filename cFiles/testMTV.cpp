/*
compile and test with the following command:
g++ -o testMVT.out -I testMVT.cpp MVTClass.hpp && ./testMVT.out
*/

#include <stdio.h>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "MVTClass.hpp"


class CSVRow
{
    public:
        std::string_view operator[](std::size_t index) const
        {
            return std::string_view(&m_line[m_data[index] + 1], m_data[index + 1] -  (m_data[index] + 1));
        }
        std::size_t size() const
        {
            return m_data.size() - 1;
        }
        void readNextRow(std::istream& str)
        {
            std::getline(str, m_line);

            m_data.clear();
            m_data.emplace_back(-1);
            std::string::size_type pos = 0;
            while((pos = m_line.find(',', pos)) != std::string::npos)
            {
                m_data.emplace_back(pos);
                ++pos;
            }
            // This checks for a trailing comma with no data after it.
            pos   = m_line.size();
            m_data.emplace_back(pos);
        }
    private:
        std::string         m_line;
        std::vector<int>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}   


class CSVIterator
{   
    public:
        typedef std::input_iterator_tag     iterator_category;
        typedef CSVRow                      value_type;
        typedef std::size_t                 difference_type;
        typedef CSVRow*                     pointer;
        typedef CSVRow&                     reference;

        CSVIterator(std::istream& str)  :m_str(str.good()?&str:nullptr) { ++(*this); }
        CSVIterator()                   :m_str(nullptr) {}

        // Pre Increment
        CSVIterator& operator++()               {if (m_str) { if (!((*m_str) >> m_row)){m_str = nullptr;}}return *this;}
        // Post increment
        CSVIterator operator++(int)             {CSVIterator    tmp(*this);++(*this);return tmp;}
        CSVRow const& operator*()   const       {return m_row;}
        CSVRow const* operator->()  const       {return &m_row;}

        bool operator==(CSVIterator const& rhs) {return ((this == &rhs) || ((this->m_str == nullptr) && (rhs.m_str == nullptr)));}
        bool operator!=(CSVIterator const& rhs) {return !((*this) == rhs);}
    private:
        std::istream*       m_str;
        CSVRow              m_row;
};


class CSVRange
{
    std::istream&   stream;
    public:
        CSVRange(std::istream& str)
            : stream(str)
        {}
        CSVIterator begin() const {return CSVIterator{stream};}
        CSVIterator end()   const {return CSVIterator{};}
};

std::vector<std::vector<std::string>> readCSVFile
(
    const char* filePath
)
{
    std::ifstream       file(filePath);

    std::vector<std::vector<std::string>> data;

    for(auto& row: CSVRange(file))
    {
        std::vector<std::string> row_data;
        for(std::size_t i = 0; i < row.size(); ++i)
        {
            row_data.push_back(std::string(row[i]));
        }
        data.push_back(row_data);
    }

    return data;
}

// function to test the MVTAnalysis class on the command line
void testMVTAnalysis(
    const char* filePath
)
{
    auto data = readCSVFile(filePath);
    // assign the data to the appropriate vectors
    std::vector<double> time;
    std::vector<double> rate;
    std::vector<double> rateErr;
    for (int i = 1; i < data.size(); i++)
    {
        time.push_back(std::stod(data[i][0]));
        rate.push_back(std::stod(data[i][1]));
        rateErr.push_back(std::stod(data[i][2]));
    }


    // convert the vectors to arrays
    double timeArray[time.size()];
    double rateArray[rate.size()];
    double rateErrArray[rateErr.size()];
    for (int i = 0; i < time.size(); i++)
    {
        timeArray[i] = time[i];
        rateArray[i] = rate[i];
        rateErrArray[i] = rateErr[i];
    }

    // print the data to the console
    for (int i = 0; i < time.size(); i++)
    {
        printf("%f\t%f\t%f\n", timeArray[i], rateArray[i], rateErrArray[i]);
    }
/*
    // test the MVTAnalysis class
    MVTAnalysis analysis(
        rateArray,
        timeArray,
        rateErrArray,
        rate.size()
    );


    // print the k set and VT set
    for (int i = 0; i < analysis.kSet.size(); i++)
    {
        std::cout << "k: " << analysis.kSet[i] << ", VT: " << analysis.VTSet[i] << std::endl;
    }
*/
}

int main()
{
    // test the MVTAnalysis class on the command line
    const char* filePath = "testData.csv";
    printf("Testing MVTAnalysis class on file: %s\n", filePath);
    testMVTAnalysis(filePath);
    return 0;
}