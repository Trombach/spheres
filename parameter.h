#ifndef PARAMETER
#define PARAMETER

#include <string>
#include <map>

template <typename T> class parameter
{
    
    private:
        std::map<std::string, T> parameters; 

    public:
        parameter() {};
        parameter(std::string _s, T _param) {parameters[_s] = _param;}

        void set(std::string _s, T _param) {parameters[_s] = _param;}

        const T get(std::string _s) {return parameters.find(_s)->second;}
};

#endif
