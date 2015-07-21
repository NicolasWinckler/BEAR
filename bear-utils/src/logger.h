/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             * 
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *  
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
/**
 * logger.h
 *
 * @since 2012-12-04
 * @author D. Klein, A. Rybalchenko
 */
// refactored from FairRoot::FairMQ project

#ifndef FAIRMQLOGGER_H_
#define FAIRMQLOGGER_H_

#include <sstream>
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <string>
#include <stdio.h>

/*
 
 #include "colormod.h" // namespace Color
#include <iostream>
using namespace std;
int main() {
    Color::Modifier red(Color::FG_RED);
    Color::Modifier def(Color::FG_DEFAULT);
    cout << "This ->" << red << "word" << def << "<- is red." << endl;
}
 */

namespace bear
{
    
    class logger
    {
      public:
        enum Level
        {
            MAXDEBUG = -1,
            DEBUG,
            INFO,
            WARNING,
            ERROR,
            NOLOG
        };

        logger();
        virtual ~logger();
        std::ostringstream& Log(int type);

        static void SetLogLevel(int logLevel)
        {
            logger::fMinLogLevel = logLevel;
        }

      private:
        std::ostringstream os;
        int fLogLevel;
        static int fMinLogLevel;
    };
    
    
    // Helper functions and macros definitions
    
    typedef unsigned long long timestamp_t;
    timestamp_t get_timestamp();
    
    #define LOG(type) bear::logger().Log(bear::logger::type)
    #define SET_LOG_LEVEL(loglevel) bear::logger::SetLogLevel(bear::logger::loglevel)
    #define SET_LOGGER_LEVEL(loglevel) bear::logger::SetLogLevel(loglevel)

    namespace color 
    {
        enum code 
        {
            FG_BLACK    = 30,
            FG_RED      = 31,
            FG_GREEN    = 32,
            FG_YELLOW   = 33,
            FG_BLUE     = 34,
            FG_MAGENTA  = 35,
            FG_CYAN     = 36,
            FG_WHITE    = 37,
            FG_DEFAULT  = 39,
            BG_RED      = 41,
            BG_GREEN    = 42,
            BG_BLUE     = 44,
            BG_DEFAULT  = 49
        };
    }
    
    template <color::code color> 
    std::string write_in(const std::string& text_in_bold)
    {
        std::ostringstream os;
        
        os << "\033[01;" << color << "m" << text_in_bold << "\033[0m";
        
        return os.str();
    }
    
    
    /*
    struct lvl 
    {
        static const logger::Level  MAXDEBUG = logger::MAXDEBUG;
        static const logger::Level  DEBUG    = logger::DEBUG;
        static const logger::Level  INFO     = logger::INFO;
        static const logger::Level  WARNING  = logger::WARNING;
        static const logger::Level  ERROR    = logger::ERROR;
        static const logger::Level  NOLOG    = logger::NOLOG;
    };
    
        
    template <logger::Level level> 
    struct log_format;
    
    template <> 
    struct log_format<logger::MAXDEBUG>
    {
        log_format(std::string& val){val=write_in<color::FG_CYAN>("MAXDEBUG");}
    };
    
    template <> 
    struct log_format<logger::DEBUG>
    {
        log_format(std::string& val){val=write_in<color::FG_BLUE>("DEBUG");}
    };

    template <> 
    struct log_format<logger::INFO>
    {
        log_format(std::string& val){val=write_in<color::FG_GREEN>("INFO");}
    };

    template <> 
    struct log_format<logger::WARNING>
    {
        log_format(std::string& val){val=write_in<color::FG_YELLOW>("WARNING");}
    };

    template <> 
    struct log_format<logger::ERROR>
    {
        log_format(std::string& val){val=write_in<color::FG_RED>("ERROR");}
    };

    template <> 
    struct log_format<logger::NOLOG>
    {
        log_format(std::string& val){val=write_in<color::FG_DEFAULT>("NOLOG");}
    };
    //*/
}



#endif /* FAIRMQLOGGER_H_ */
