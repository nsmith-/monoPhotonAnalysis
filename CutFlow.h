#ifndef CUTFLOW_H
#define CUTFLOW_H

#include <vector>
#include <map>
#include <ostream>

class CutFlow
{
  public:
    CutFlow();
    virtual ~CutFlow();

    void increment(std::string label) { if ( counts[label]++ == 1 ) order.push_back(label); };
    friend std::ostream& operator<<(std::ostream& os, CutFlow& flow)
    {
      for(const auto label : flow.order)
      {
        os << std::setw(30) << label << " : " << flow.counts[label] << " events passed." << std::endl;
      }
      return os;
    };

  private:
    std::map<std::string, long> counts;
    std::vector<std::string> order;
};

#endif
