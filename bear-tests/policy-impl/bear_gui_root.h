/* 
 * File:   bear_gui_root.h
 * Author: winckler
 *
 * Created on August 10, 2015, 6:36 PM
 */

#ifndef BEAR_GUI_ROOT_H
#define	BEAR_GUI_ROOT_H

#include <map>
#include <vector>
#include <string>
#include <memory>
#include <iostream>     
#include <sstream>      // std::stringstream
#include <type_traits>
#include <limits>       // std::numeric_limits

#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"

#include "logger.h"
#include "def.h"

namespace bear
{
    class bear_gui_root
    {
        
    public:
        
        enum method {kDiagonalization,kRungeKutta};
        
        bear_gui_root() :   fCanvas(nullptr), 
                            fLegend(nullptr),   
                            fFunctions(),   
                            fXmin(0.),  
                            fXmax(20.),     
                            fYmin(0.),  
                            fYmax(1.1),
                            fNpoint(1000),
                            fMaxop(100000),
                            fMaxpar(100000),
                            fMaxconst(100000),
                            fFunctions_derivative()
        {
        }
        virtual ~bear_gui_root()
        {
            //for(auto& p : fFunctions)
                //if(p.second)
                   // delete p.second;
            
            if(fLegend)
                delete fLegend;
            
            if(fCanvas)
                delete fCanvas;
            
            
        }
        
        int init(const variables_map& vm,const variables_map& vm2)
        {
            
            fXmin=vm.at("thickness.minimum").template as<double>();
            fXmax=vm.at("thickness.maximum").template as<double>();
            fNpoint=vm.at("thickness.point.number").template as<std::size_t>();
            
            fYmin=vm.at("fraction.minimum").template as<double>();
            fYmax=vm.at("fraction.maximum").template as<double>();
            
            std::string proj_symbol=vm.at("projectile.symbol").template as<std::string>();
            std::string proj_energy=vm.at("projectile.energy").template as<std::string>();
            std::string target_symbol=vm.at("target.symbol").template as<std::string>();
            int tmass=(int)vm.at("target.mass.number").template as<double>();
            std::string target_mass=std::to_string(tmass);
            std::string target_pressure=vm.at("target.pressure").template as<std::string>();
            std::string X_unit=vm.at("thickness.unit").template as<std::string>();
            
            std::stringstream ss;
            
            ss<<proj_symbol
                    <<" projectile at "
                    <<proj_energy
                    <<" on ^{"
                    <<target_mass<<"}"
                    <<target_symbol
                    <<" target"
                    ;
            ss<<" with "<< target_pressure <<" pressure.";
            
            fTitle=ss.str();
            fXTitle="thickness (";
            if(X_unit!="mug/cm2")
                fXTitle+=X_unit;
            else
                fXTitle+="#mug/cm2";
            fXTitle+=")";
            
            
            fYTitle="Fractions";
            
            fLegend = new TLegend(0.4,0.7,0.9,0.9);
            fLegend->SetNColumns(4);

            LOG(DEBUG)<<"init(variable_map)";
            fMaxop=vm2.at("formula-maximum-operator").template as<int>();
            fMaxpar=vm2.at("formula-maximum-parameter").template as<int>();
            fMaxconst=vm2.at("formula-maximum-constant").template as<int>();
            return 0;
        }

        
        int init_summary(std::shared_ptr<bear_summary> const& summary) 
        {
            fSummary = summary;
            return 0;
        }
        
        
        int compute_equilibrium_distance()
        {
            double epsilon=0.0001;
            double x=0;
            for(auto& p : fFunctions)
            {
                double val=p.second->Derivative(x,nullptr,epsilon);
            }
                
            return 0;
        }
        
        int print_table()
        {
            
            std::ostringstream os_title;
            os_title<< std::setw(16)
                    //<<std::left
                    << bstream_centered("X") << "    ";
            for(auto& p : fFunctions)
            {
                os_title<< std::setw(16)
                        //<<std::left
                        << bstream_centered("F"+std::to_string( p.first + 1)) << "    ";
            }
            os_title<<std::setw(16)
                    //<<std::left
                    <<bstream_centered("Sum");
            LOG(RESULTS)<<os_title.str();
            double N=(double)fNpoint;
            double step =(fXmax-fXmin)/N;
            double sum=0;
            for(int i(0);i<fNpoint;i++)
            {
                sum=0;
                std::ostringstream os_eval;
                double idoub = (double)i;
                double x=idoub*step+fXmin;
                os_eval<< std::setw(16)
                        //<<std::left
                        //<< std::setprecision(10) 
                        << std::scientific
                        <<bstream_centered(to_string_scientific(x))<<"    ";
                for(auto& p : fFunctions)
                {
                    
                    os_eval << std::setw(16)
                            //<<std::left
                            //<< std::setprecision(12) 
                            << bstream_centered(to_string_scientific(p.second->Eval(x))) << "    ";
                    sum+=p.second->Eval(x);
                }
                os_eval << std::setw(16)
                        //<< std::left 
                        << bstream_centered(to_string_scientific(sum));
                LOG(RESULTS)<<os_eval.str();
            }
            return 0;
        }
        
        
        
        // init functions/histos
        int init(const std::map<std::size_t, std::string>& input_functions)
        {
            //Int_t maxop=std::numeric_limits<Int_t>::max();
            //Int_t maxpar=std::numeric_limits<Int_t>::max();
            //Int_t maxconst=std::numeric_limits<Int_t>::max();
            
            
            fMethod=kDiagonalization;
            for(const auto& p : input_functions)
            {
                std::string name = "F" + std::to_string(fSummary->F_index_map.at(p.first));
                fFunctions[p.first] = std::make_shared<TF1>(name.c_str(), p.second.c_str(), fXmin, fXmax);
                fFunctions.at(p.first)->SetNpx(fNpoint);
                fFunctions.at(p.first)->SetLineColor(p.first+1);
                fFunctions.at(p.first)->SetMaxima(fMaxop,fMaxpar,fMaxconst);
                fLegend->AddEntry(fFunctions[p.first].get(), name.c_str());
            }
            return 0;
        }
        
        int init(std::map<std::size_t, std::shared_ptr<TH1D> >& input_functions, bool plot=false)
        {
            fMethod=kRungeKutta;
            fTitle+=" (Runge Kutta method)";
            fHistograms=input_functions;
            for(auto& p : fHistograms)
            {
                std::string name = "F" + std::to_string(p.first+1);
                p.second->SetLineColor(p.first+1);
                p.second->SetLineWidth(2);
                p.second->SetStats(kFALSE);
                fLegend->AddEntry(p.second.get(), name.c_str());
            }
            //temporary hack
            if(plot)
                return draw(fHistograms);
            return 0;
        }
        
        // plot and draw functions
        int plot()
        {
            if(fMethod==kDiagonalization)
                return draw(fFunctions);
            
            //temporary hack
            if(fMethod==kRungeKutta)
                return init(fHistograms,true);
            
            
            return 0;
        }
        
        
        
        template <typename T>
        int draw(std::map<std::size_t,T>& container_map)
        {
            LOG(DEBUG)<<"GUI start";
            
            if(fMethod==kDiagonalization)
                fCanvas = new TCanvas("c1Dia","Solutions - Diagonalization",800,600);
            
            
            if(fMethod==kRungeKutta)
                fCanvas = new TCanvas("c1RK","Solutions - Runge Kutta",800,600);
            
            
            if(fMethod!=kRungeKutta && fMethod!=kDiagonalization)
                throw std::runtime_error("Unrecognized method to solve the equations");
            
            fCanvas->SetLogx();
            for(auto& p : container_map)
            {
                
                if(p.first!=fSummary->max_fraction_index)
                    p.second->Draw("SAME");
                else
                {
                    p.second->SetTitle(fTitle.c_str());
                    p.second->GetYaxis()->SetRangeUser(fYmin,fYmax);
                    
                    
                    p.second->GetXaxis()->CenterTitle();
                    p.second->GetYaxis()->CenterTitle();
                    
                    p.second->GetXaxis()->SetTitleOffset(1.2);
                    p.second->GetYaxis()->SetTitleOffset(1.2);
                    
                    p.second->GetXaxis()->SetTitle(fXTitle.c_str());
                    p.second->GetYaxis()->SetTitle(fYTitle.c_str());
                    
                    p.second->Draw();
                }
            }

            //fCanvas->SetLogx();
            fLegend->Draw();
            return 0;
        }
        
    private:
        TCanvas* fCanvas;
        TLegend* fLegend;
        //TF1 *fa1;
        std::map<std::size_t, std::string> fInput;
        //std::map<std::size_t, TF1*> fFunctions;
        std::map<std::size_t, std::shared_ptr<TF1> > fFunctions;
        std::map<std::size_t, std::shared_ptr<TF1> > fFunctions_derivative;
        double fXmin;
        double fXmax;
        double fYmin;
        double fYmax;
        std::string fTitle;
        std::string fXTitle;
        std::string fYTitle;
        std::size_t fNpoint;
        enum method fMethod;
        std::shared_ptr<bear_summary> fSummary;
        std::map<std::size_t, std::shared_ptr<TH1D> > fHistograms;
        Int_t fMaxop;
        Int_t fMaxpar;
        Int_t fMaxconst;
    };
}


#endif	/* BEAR_GUI_ROOT_H */

