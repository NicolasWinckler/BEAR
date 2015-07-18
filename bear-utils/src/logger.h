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

namespace bear
{
    class logger
    {
      public:
        enum Level
        {
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

    typedef unsigned long long timestamp_t;

    timestamp_t get_timestamp();

    
}

#define LOG(type) bear::logger().Log(bear::logger::type)
#define SET_LOG_LEVEL(loglevel) bear::logger::SetLogLevel(bear::logger::loglevel)
#define SET_LOGGER_LEVEL(loglevel) bear::logger::SetLogLevel(loglevel)

#endif /* FAIRMQLOGGER_H_ */
