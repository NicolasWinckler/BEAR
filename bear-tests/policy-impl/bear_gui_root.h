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
#include "TFormula.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"


#include "logger.h"
#include "def.h"
#include "handle_root_signal.h"

namespace bear
{
    class bear_gui_root
    {
        
    public:
        
        enum method {kDiagonalization,kRungeKutta};
        
        bear_gui_root() :   fCanvas_non_equilib(nullptr), 
                            fCanvas_equilib(nullptr),
                            fLegend(nullptr),
                            fLegend_eq(nullptr),
                            fEquilibrium_solutions(nullptr),
                            fFunctions(),
                            fXmin(0.),  
                            fXmax(20.),     
                            fYmin(0.),  
                            fYmax(1.1),
                            fNpoint(1000),
                            fMaxop(100000),
                            fMaxpar(100000),
                            fMaxconst(100000),
                            fFunctions_derivative(),
                            fSingal_handler(),
                            fOut_fig_filename(),
                            fOut_fig_e_filename(),
                            fCharge(nullptr),
                            fFraction(nullptr),
                            fSave_ne(false),
                            fSave_e(false)
        {
        }
        virtual ~bear_gui_root()
        {
            fCanvas_non_equilib.reset();
            fCanvas_equilib.reset();
            fLegend.reset();
            fLegend_eq.reset();
            fEquilibrium_solutions.reset();
            delete[] fCharge;
            delete[] fFraction;
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
            
            fLegend = std::make_shared<TLegend>(0.4,0.7,0.9,0.9);
            fLegend->SetNColumns(4);
            fLegend_eq = std::make_shared<TLegend>(0.1,0.7,0.4,0.9);

            LOG(DEBUG)<<"init(variable_map)";
            fMaxop=vm2.at("formula-maximum-operator").template as<int>();
            fMaxpar=vm2.at("formula-maximum-parameter").template as<int>();
            fMaxconst=vm2.at("formula-maximum-constant").template as<int>();
            
            fs::path input=vm2["input-file"].template as<fs::path>();
            std::string filename=input.stem().string();
            std::string outputdir=vm2["output-directory"].template as<fs::path>().string();
            fOut_fig_e_filename=outputdir;
            fOut_fig_filename=outputdir;

            fOut_fig_filename+="/Bear-results-figure-ne-";
            fOut_fig_filename+=filename;// pdf/root extension added later

            fOut_fig_e_filename+="/Bear-results-figure-e-";
            fOut_fig_e_filename+=filename;// pdf/root extension added later

            //output+=".root";
            //LOG(INFO)<<"save figure to : "<<output;
            
            fSave_ne=vm2["save-fig-ne"].template as<bool>();
            fSave_e=vm2["save-fig-e"].template as<bool>();
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
                        << bstream_centered("F"+std::to_string( fSummary->F_index_map.at(p.first) )) << "    ";
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
            //fCanvas_non_equilib = std::make_shared<TCanvas>("c1Dia","Solutions - Diagonalization",800,600);
            
            int color=1;
            for(const auto& p : input_functions)
            {
                std::string name = "F" + std::to_string(fSummary->F_index_map.at(p.first));
                fFunctions[p.first] = std::make_shared<TF1>(name.c_str(), p.second.c_str(), fXmin, fXmax);
                fFunctions.at(p.first)->SetNpx(fNpoint);
                int color =p.first;
                if(color == 0 || color == 10)
                    color++;

                if(color==50)
                    color=1;

                fFunctions.at(p.first)->SetLineColor(color);
                #ifdef __CINT__
                fFunctions.at(p.first)->SetMaxima(fMaxop,fMaxpar,fMaxconst); // deprecated in root 6
                #endif
                fLegend->AddEntry(fFunctions[p.first].get(), name.c_str());
            }
            return 0;
        }
        
        int init(std::map<std::size_t, std::shared_ptr<TH1D> >& input_functions, bool plot=false)
        {
            fMethod=kRungeKutta;
            fTitle+=" (Runge Kutta method)";
            fHistograms=input_functions;
            fCanvas_non_equilib = std::make_shared<TCanvas>("c1RK","Solutions - Runge Kutta",800,600);
            
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
        
        int save_fig(const std::string& filename)
        {
            fCanvas_non_equilib->SaveAs(filename.c_str());
            return 0;
        }
        
        int save_fig_equilibrium(const std::string& filename)
        {
            fCanvas_equilib->SaveAs(filename.c_str());
            return 0;
        }
        
        template <typename T>
        int draw(std::map<std::size_t,T>& container_map)
        {
            LOG(DEBUG)<<"GUI start";
            
            fCanvas_non_equilib = std::make_shared<TCanvas>("c1Dia","Solutions - Diagonalization",800,600);
            
            
            if(fMethod!=kRungeKutta && fMethod!=kDiagonalization)
                throw std::runtime_error("Unrecognized method to solve the equations");
            
            fSingal_handler.set_canvas(fCanvas_non_equilib.get());

            fCanvas_non_equilib->SetLogx();

            // need two loops for draw and draw same
            for(auto& p : container_map)
            {
                if(p.first==fSummary->max_fraction_index)
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

            for(auto& p : container_map)
            {
                if(p.first!=fSummary->max_fraction_index)
                    p.second->Draw("SAME");
            }


            //fCanvas_non_equilib->SetLogx();
            fLegend->Draw();
            
            
            if(fSave_ne)
            {
                LOG(INFO)<<" ";
                LOG(STATE)<<"saving to figure ...";
                LOG(INFO)<<"- saving figure of non-equilibrium solutions to : ";//<<fOut_fig_filename;
                save_fig(fOut_fig_filename+".pdf");
                save_fig(fOut_fig_filename+".root");
                LOG(INFO)<<" ";
            }


            bool plot_eq=true;
            if(plot_eq)
            {
                std::size_t dim=fSummary->equilibrium_solutions.size(); 
                const Int_t max_dim=15;
                Double_t x[max_dim], y[max_dim];
                fCharge = new double[dim];
                fFraction = new double[dim];

                
                if(dim>max_dim)
                {
                    LOG(ERROR)<< "The dimension of the system is limited to N=200";
                    return 1;
                }

                if(dim==0)
                {
                    LOG(ERROR)<< "The dimension of the system is zero, check input file";
                    return 1;
                }

                std::size_t i(0);
                double charge_start;
                double charge_end;
                double mean=0.;
                double sum=0.;
                for(const auto& p : fSummary->equilibrium_solutions)
                {
                    double Fi=p.second;
                    double qi=fSummary->F_index_map.at(p.first);


                    fCharge[i]=qi;
                    fFraction[i]=Fi;
                    mean+=qi*Fi;
                    i++;
                    //mean+=qi*Fi;
                    sum+=Fi;
                    //LOG(ERROR)  << "F"
                    //              << std::to_string(fSummary->F_index_map.at(p.first))
                    //              << " = "
                    //              << Fi;
                }

                double variance=0.;
                for(const auto& p : fSummary->equilibrium_solutions)
                {
                    double Fi=p.second;
                    double qi=fSummary->F_index_map.at(p.first);
                    double temp=(qi-mean);
                    int N=static_cast<int>(dim);
                    variance+=temp*temp/N;
                }

                double sigma=std::sqrt(variance);


                charge_start=fCharge[0]-1;
                charge_end=fCharge[dim-1]+1;

                fCanvas_equilib = std::make_shared<TCanvas>("c1equi","Solutions at equilibrium",800,600);
                fSingal_handler.set_canvas2(fCanvas_equilib.get());
                fEquilibrium_solutions = std::make_shared<TGraph>(dim,fCharge,fFraction);

                fEquilibrium_solutions->SetLineColor(4);
                fEquilibrium_solutions->SetLineWidth(2);
                fEquilibrium_solutions->SetMarkerColor(2);
                fEquilibrium_solutions->SetMarkerSize(1.5);
                fEquilibrium_solutions->SetMarkerStyle(20);
                std::string graph_title=fTitle+" (equilibrium)";
                fEquilibrium_solutions->SetTitle(graph_title.c_str());


                fEquilibrium_solutions->GetXaxis()->CenterTitle();
                fEquilibrium_solutions->GetYaxis()->CenterTitle();
                
                fEquilibrium_solutions->GetXaxis()->SetTitleOffset(1.2);
                fEquilibrium_solutions->GetYaxis()->SetTitleOffset(1.2);
                fEquilibrium_solutions->GetXaxis()->SetTitle("Charge q");
                fEquilibrium_solutions->GetXaxis()->SetRangeUser(charge_start,charge_end);

                fEquilibrium_solutions->GetYaxis()->SetTitle("Fractions");

                fEquilibrium_solutions->Draw("ALP");
                std::string graph_legend("#splitline{Equilibrium charge}{state distribution}");
                
                fLegend_eq->AddEntry((TObject*)0,"","");
                fLegend_eq->AddEntry(fEquilibrium_solutions.get(),graph_legend.c_str(),"lp");
                
                fLegend_eq->AddEntry((TObject*)0,"","");
                graph_legend="<q> = "+std::to_string(mean);
                fLegend_eq->AddEntry((TObject*)0,graph_legend.c_str(),"");

                fLegend_eq->AddEntry((TObject*)0,"","");
                graph_legend="#sigma = "+std::to_string(sigma);
                fLegend_eq->AddEntry((TObject*)0,graph_legend.c_str(),"");
                fLegend_eq->AddEntry((TObject*)0,"","");
                
                for(const auto& p : fSummary->equilibrium_solutions)
                {
                    double Fi=p.second;
                    double qi=fSummary->F_index_map.at(p.first);
                    //mean+=qi*Fi;
                    //sum+=Fi;
                    LOG(DEBUG)    << "mean = "<< mean << " F"
                                  << std::to_string(fSummary->F_index_map.at(p.first))
                                  << " = "
                                  << Fi;
                }

                fLegend_eq->Draw();

                if(fSave_e)
                {
                    
                    if(!fSave_ne)
                    {
                        LOG(INFO)<<" ";
                        LOG(STATE)<<"saving to figure ...";
                    }

                    LOG(INFO)<<"- saving figure of equilibrium solutions to : ";//<<fOut_fig_e_filename;
                    save_fig(fOut_fig_e_filename+".pdf");
                    save_fig(fOut_fig_e_filename+".root");
                    LOG(INFO)<<" ";
                }
            }



            return 0;
        }
        
    private:
        std::shared_ptr<TCanvas> fCanvas_non_equilib;
        std::shared_ptr<TCanvas> fCanvas_equilib;
        std::shared_ptr<TLegend> fLegend;
        std::shared_ptr<TLegend> fLegend_eq;
        std::shared_ptr<TGraph>  fEquilibrium_solutions;
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
        
        handle_root_signal fSingal_handler;
        std::string fOut_fig_filename;
        std::string fOut_fig_e_filename;
        
        double* fCharge;     // for equilibrium solution
        double* fFraction;   // for equilibrium solution 
        bool fSave_ne;
        bool fSave_e;
        // 
    };
}


#endif	/* BEAR_GUI_ROOT_H */

