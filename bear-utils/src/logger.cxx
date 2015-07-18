/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             * 
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *  
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
/**
 * logger.cxx
 *
 * @since 2012-12-04
 * @author D. Klein, A. Rybalchenko
 */

// refactored from FairRoot::FairMQ project


#include "logger.h"

namespace bear
{
    int logger::fMinLogLevel = logger::DEBUG;

    logger::logger()
        : os()
        , fLogLevel(DEBUG)
    {
    }

    logger::~logger()
    {
        if (fLogLevel >= logger::fMinLogLevel && fLogLevel < logger::NOLOG)
        {
            std::cout << os.str() << std::endl;
        }
    }

    std::ostringstream& logger::Log(int type)
    {
        std::string type_str;
        fLogLevel = type;
        switch (type)
        {
            
            case DEBUG :
                type_str = "\033[01;34mDEBUG\033[0m";
                break;

            case INFO :
                type_str = "\033[01;32mINFO\033[0m";
                break;

            case WARNING :
                type_str = "\033[01;33mWARN\033[0m";
                break;

            case ERROR :
                type_str = "\033[01;31mERROR\033[0m";
                break;    

            //case STATE :
            //    type_str = "\033[01;35mSTATE\033[0m";
            //    break;

            case NOLOG :
                type_str = "\033[01;31mNOLOG\033[0m";
                break;

            default:
                break;
        }

        timestamp_t tm = get_timestamp();
        timestamp_t ms = tm / 1000.0L;
        timestamp_t s = ms / 1000.0L;
        time_t t = s;
        // size_t fractional_seconds = ms % 1000;
        char mbstr[100];
        strftime(mbstr, 100, "%H:%M:%S", localtime(&t));

        os << "[\033[01;36m" << mbstr << "\033[0m]"
           << "[" << type_str << "]"
           << " ";

        return os;
    }

    timestamp_t get_timestamp()
    {
        struct timeval now;
        gettimeofday(&now, NULL);
        return now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
    }

}// end namespace