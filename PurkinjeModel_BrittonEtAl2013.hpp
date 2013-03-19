/*
 * PurkinjeModel_BrittonEtAl2013.hpp
 *
 *  Created on: 4 Mar 2013
 *      Author: olibri
 */

/*

Copyright (C) University of Oxford, 2005-2013

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef PURKINJEMODEL_BRITTONETAL2013_HPP_
#define PURKINJEMODEL_BRITTONETAL2013_HPP_

#include <iostream>
#include <fstream>
#include <sstream>

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractModelPopCvodeSystem.hpp"
#include "VectorHelperFunctions.hpp"
#include "OdeSystemInformation.hpp"



class PurkinjeModel_BrittonEtAl2013 : public AbstractModelPopCvodeSystem
{

public:
	PurkinjeModel_BrittonEtAl2013() :
		AbstractModelPopCvodeSystem(29)
	{
		mpSystemInfo = OdeSystemInformation<PurkinjeModel_BrittonEtAl2013>::Instance();
		Init();
	}

	//Constructor with initial conditions file specified
	PurkinjeModel_BrittonEtAl2013(std::string initialConditionsFileName) :
		AbstractModelPopCvodeSystem(29)
	{
		mpSystemInfo = OdeSystemInformation<PurkinjeModel_BrittonEtAl2013>::Instance();
		LoadInitialConditions(initialConditionsFileName);
		Init();
	}

	void EvaluateYDerivatives(double time, const N_Vector N_Y, N_Vector N_dY)
	{
		std::vector<double> Y;
		std::vector<double> dY;

		//Convert Y from N_Vector to std::vectors
		double tempY;
		for (int i=0; i<29; i++)
		{
			tempY = GetVectorComponent(N_Y,i);
			Y.push_back(tempY);
			dY.push_back(0.0);
		}

		// ---------
		// CONSTANTS
		// ---------

		// Constants from Parameter Set

		double G_max_na_fast = mParameterSet[0];
		double G_max_na_late = mParameterSet[1];
		double G_max_Ltype = mParameterSet[2];
		double G_max_kr = mParameterSet[3];
		double G_max_ks = mParameterSet[4];
		double G_max_k1 = mParameterSet[5];
		double G_max_to_fast = mParameterSet[6];
		double G_max_to_sustained = mParameterSet[7];
		double tau_x_Ltype = mParameterSet[8];
		double tau_x_na_fast = mParameterSet[9];
		double tau_x_na_late = mParameterSet[10];
		double g_nak = mParameterSet[11];

		// Constants from Corrias Model

		double F=96485.3415; // Faraday Constant (C/M)
	    double R=8314.472; // Gas Constant (J/(K mM))
	    double T=310.0; // Temperature (K)
	    double Ca_o = 2.0;   // millimolar (in Environment)
		double K_o = 4.0;// millimolar (in Environment)
		// K_o changed from 5.4 mM to match experimental conditions
		double Na_o = 140.0;   // millimolar (in Environment)
		double G_max_cab = 0.0001;   // conductance_units (in i_cab)
		double G_max_Ttype = 0.9;   // conductance_units (in i_cat)
		double G_f_k = 0.188709677;   // conductance_units (in i_f_k)
		double G_f_na = 0.045290323;   // conductance_units (in i_f_na)
		double G_max_kb = 0.01;   // conductance_units (in i_kb)
		double G_max_nab = 0.01;   // conductance_units (in i_nab)
		double Hpmca = 1.5;   // dimensionless (in i_pmca)
		double Kpmca = 0.0001;   // millimolar (in i_pmca)
		double PMCA_max = 5.0;   // current_units (in i_pmca)

		double Cm = 69.0;   // capacitance_units (in membrane)
		double bulk_fraction = 0.6;   // dimensionless (in membrane)
		double cell_volume = 13266.5;   // volume_units (in membrane)
		double periphery_fraction = 0.17;   // dimensionless (in membrane)
		// Changed to 0.17 to be consistent with PRd changes

		double tau_x_Ttype = 1.0;   // time_units (in x_Ttype)
		double tau_x_to_fast = 5.0;   // time_units (in x_to_fast)
		double tau_y_to_fast = 350.0;   // time_units (in y_to_fast)

		// Stimulus related constants

		double stim_amplitude = -2300.0;//-4320.0;   // current_units (in membrane)
		double stim_duration = 2.0;   // time_units (in membrane)
		double stim_period = mParameterSet.back();   // time_units (in membrane)
		double stim_start = 0.0;   // time_units (in membrane)
		double stim_end = 1000000000.0;   // time_units (in membrane)

		// Constants from PRd model

		double v_max=2.925;

		double ksat=0.27, eta=0.35, KmNai=12.3, KmNao=87.5, KmCai=0.0036, KmCao=1.3, KmCa_act=0.000125;

		double k0=96000.0, k0a=9.6, k1=150000.0, k1a=16.5, k2=1800.0, k2a=0.21, t_ip3r=3.7, IP3=0.001;

		double d_K_plb=0.00017, delta_J=0.75, Km_camk=0.15, Jserca=0.0026, Jserca_s=0.0002, Km_serca=0.00028, NSR_avg=15.0;

		double t_tr=120.0;

		double t_diff=0.2, t_gap=12.0;

		double alpha_camk=0.05, beta_camk=0.00068, CAMK0=0.05, Km_cam=0.0015;

		double BSRbar = 0.02, KmBSR = 0.00087;
		double BSLbar = 0.48, KmBSL = 0.0087;

		double V_myo=0.60, V_nsr=0.04, V_jsr=0.002, V_csr=0.008, V_pcs=0.02, V_ssl=0.15;


		// ------------------
		// COMPUTED VARIABLES
		// ------------------

		double RToF = R*T/F;
		double diffusable_volume = (bulk_fraction+periphery_fraction)*cell_volume;

		double E_ca = 0.5*RToF*log(Ca_o/Y[1]);
		double i_cab = G_max_cab*(Y[5]-E_ca);
		double E_Ca_1 = 0.5*RToF*log(Ca_o/Y[24]); //Change to 24 - in PCS not SSL!
		double i_cal = G_max_Ltype*Y[6]*Y[12]*Y[14]*(Y[5]-E_Ca_1);
		double E_Ca_2 = 0.5*RToF*log(Ca_o/Y[1]);
		double i_cat = G_max_Ttype*Y[7]*Y[13]*(Y[5]-E_Ca_2);
		double E_k_1 = RToF*log(K_o/Y[3]);
		double i_f_k = G_f_k*Y[15]*(Y[5]-E_k_1);
		double E_na_1 = RToF*log(Na_o/Y[4]);
		double i_f_na = G_f_na*Y[16]*(Y[5]-E_na_1);
		double E_k_2 = RToF*log(K_o/Y[3]);
		double x_k1 = 1.0/(1.0+exp((92.0+Y[5])/10.0));
		double i_k1 = G_max_k1*pow((K_o/5.4),0.8)*x_k1*(Y[5]-E_k_2);
		double E_k_3 = RToF*log(K_o/Y[3]);
		double i_kb = G_max_kb*(Y[5]-E_k_3);
		double E_k_4 = RToF*log(K_o/Y[3]);
		double x_kr = 1.0/(1.0+exp((33.0+Y[5])/22.4));
		double i_kr = G_max_kr*(K_o/5.4)*x_kr*Y[17]*(Y[5]-E_k_4);
		double E_k_5 = RToF*log(K_o/Y[3]);
		double i_ks = G_max_ks*Y[8]*Y[18]*(Y[5]-E_k_5);
		double E_na_2 = RToF*log(Na_o/Y[4]);
		double i_na_fast = G_max_na_fast*Y[9]*Y[19]*(Y[5]-E_na_2);
		double E_na_3 = RToF*log(Na_o/Y[4]);
		double i_na_late = G_max_na_late*Y[10]*Y[20]*(Y[5]-E_na_3);
		double E_na_4 = RToF*log(Na_o/Y[4]);
		double i_nab = G_max_nab*(Y[5]-E_na_4);
		double x_nak = 1.0/(1.0+exp((Y[5]+80.0)/-45.0));
		double y_nak = 1.0/(1.0+exp((Y[5]+0.0)/125.0));
		double i_nak = g_nak*x_nak*y_nak*1.0/(1.0+pow((1.9/K_o),1.45))*1.0/(1.0+(31.98/Y[4]));
		double i_pmca = PMCA_max*1.0/(1.0+pow((Kpmca/Y[1]),Hpmca));
		double E_k_6 = RToF*log(K_o/Y[3]);
		double i_to_fast = G_max_to_fast*Y[11]*Y[21]*(Y[5]-E_k_6);
		double E_k_7 = RToF*log(K_o/Y[3]);
		double x_to_sustained = 1.0/(1.0+exp((5.0-Y[5])/17.0));
		double i_to_sustained = G_max_to_sustained*x_to_sustained*(Y[5]-E_k_7);

		// Stimulus

		double i_stim;
		if ((time >= stim_start) && (time <= stim_end) && (time-stim_start-floor((time-stim_start)/stim_period)*stim_period <= stim_duration))
		{
			i_stim = stim_amplitude;
		}
		else
		{
			i_stim = 0.0;
		}

		// Computed Variables for calcium sub-system (from PRd model)

		//NACA

		double Allo_ssl = 1.0/(1.0+pow(KmCa_act/(1.5*Y[1]),2.0));
		double Allo_pcs = 1.0/(1.0+pow(KmCa_act/(1.5*Y[24]),2.0));

		double deltaE_ssl_num = (v_max*(pow(Y[4],3.0)*Ca_o*exp(eta*Y[5]*F/(R*T))-pow(Na_o,3.0)*1.5*Y[1]*exp((eta-1)*Y[5]*F/(R*T))));
		double mult = (1.0+ksat*exp((eta-1.0)*Y[5]*F/(R*T)));
		double sum1 = (KmCao*(pow(Y[4],3.0))+(pow(KmNao,3.0))*1.5*Y[1]+(pow(KmNai,3))*Ca_o*(1.0+(1.5*Y[1])/KmCai));
		double sum2 = KmCai*(pow(Na_o,3.0))*(1.0+(pow(Y[4],3.0)/pow(KmNai,3.0)))+(pow(Y[4],3.0))*Ca_o+pow(Na_o,3.0)*1.5*Y[1];

		double deltaE_ssl = deltaE_ssl_num/(mult*(sum1+sum2));

		double deltaE_pcs = (v_max*(pow(Y[4],3.0)*Ca_o*exp(eta*Y[5]*F/(R*T))-pow(Na_o,3.0)*1.5*Y[24]*exp((eta-1.0)*Y[5]*F/(R*T))))
							/ ((1.0+ksat*exp((eta-1.0)*Y[5]*F/(R*T)))*(KmCao*(pow(Y[4],3.0))+pow(KmNao,3.0)*1.5*Y[24]+pow(KmNai,3.0)*Ca_o*(1.0+(1.5*Y[24])/KmCai)+KmCai*pow(Na_o,3.0)*(1.0+(pow(Y[4],3.0)/pow(KmNai,3.0)))+pow(Y[4],3.0)*Ca_o+pow(Na_o,3.0)*1.5*Y[24]));

		double I_naca_ssl = Allo_ssl*deltaE_ssl*Cm;
		double I_naca_pcs = Allo_pcs*deltaE_pcs*Cm;
		double i_naca = (0.8*I_naca_ssl + 0.2*I_naca_pcs);

		//CAMK

		double CAMK_bound = CAMK0*(1.0-Y[25])/(1.0+(Km_cam/Y[24]));
		double CAMK_act = CAMK_bound + Y[25];

		//SERCA

		double d_Km_plb = d_K_plb*CAMK_act/(Km_camk+CAMK_act);
		double deltaJsC = delta_J*CAMK_act/(Km_camk + CAMK_act);

		double J_serca = Jserca*(1.0+deltaJsC)/(1.0+((Km_serca-d_Km_plb)/Y[0])) - 0.0042*Y[2]/NSR_avg;
		double J_serca_s = Jserca_s*(1.0+deltaJsC)/(1.0+((Km_serca-d_Km_plb)/Y[1])) - 0.00105*Y[2]/NSR_avg;

		//IP3RC

		double J_ip3r = 10.92*(Y[22] - Y[24])*(t_ip3r*IP3*Y[24]*(1.0-Y[28])/((1.0+(IP3*k0)/k0a)*(1.0+Y[24]*k1/k1a)));

		//Translocation Fluxes

		double J_tr_j = (Y[2] - Y[22])/t_tr;
		double J_tr_c = (Y[2] - Y[23])/t_tr;

		//Diffusion in cytoplasm

		double J_diff = (Y[24] - Y[1])/t_diff;
		double J_gap = (Y[1] - Y[0])/t_gap;

		//RYR3

		double Rel_ryr3 = -i_cal*1000.0/(F*V_pcs*cell_volume) + (Y[27] + J_ip3r)*(V_jsr/V_pcs) - J_diff;
		double t_ryr3 = 2.0*(1.0+(1.0/(1.0+pow((0.28/CAMK_act),8.0))))/(1.0+(0.0123/Y[22]));
		if(t_ryr3 < 0.01)
		{
			t_ryr3 = 0.01;
		}

		double Ryr3_temp = 15.0*Rel_ryr3*(1.0+(1.0/(1.0+pow((0.28/CAMK_act),8.0))))/(1.0+pow((1.0/Y[22]),8.0));
		double Ryr3_inf;
		if(Rel_ryr3 > 0.0)
		{
			Ryr3_inf=Ryr3_temp;
		}
		else
		{
			Ryr3_inf=0.0;

		}

		//RYR2

		double Rel_ryr2 = -J_serca*V_nsr/V_myo + J_gap*V_ssl/V_myo + Y[26]*V_csr/V_myo;
		double t_ryr2 = 6.0*(1.0+(1.0/(1.0+pow((0.28/CAMK_act),8.0))))/(1.0+(0.0123/Y[23]));
		double Ryr2_temp = 91.0*Rel_ryr2*(1.0+(1.0/(1.0+pow((0.28/CAMK_act),8.0))))/(1.0+pow((1.0/Y[23]),8.0));

		if(t_ryr2 < 0.01)
		{
			t_ryr2 = 0.01;
		}

		double Ryr2_inf;
		if(Rel_ryr2 > 0.0)
		{
			Ryr2_inf=Ryr2_temp;
		}
		else
		{
			Ryr2_inf=0.0;
		}

		double Beta_pcs = 1.0/(1.0 + ((BSRbar*KmBSR)/pow((Y[24] + KmBSR),2.0)) + ((BSLbar*KmBSL)/pow((Y[24] + KmBSL),2.0)));


		// --- Derivatives ---

		// Vm - Membrane Potential
		dY[5] = (-1.0*(1.0/Cm)*(i_k1+i_to_fast+i_to_sustained+i_kr+i_ks+i_kb+i_nak+i_cal+i_cat+i_na_fast+i_na_late+i_pmca+i_cab+i_f_na+i_f_k+i_nab+i_stim+i_naca));

		// - Calcium Sub-system -

		// Ca_Bulk
		dY[0] = J_gap*(V_ssl/V_myo) - J_serca*(V_nsr/V_myo) + Y[26]*(V_csr/V_myo);

		// Ca_SSL
		dY[1] = (2.0*I_naca_ssl - i_cat - i_cab - i_pmca)*1000.0/(F*diffusable_volume) + J_diff*(V_pcs/V_ssl) - J_serca_s*(V_nsr/V_ssl) - J_gap;

		// Ca_NSR
		dY[2] = J_serca + J_serca_s - J_tr_c*(V_csr/V_nsr) - J_tr_j*(V_jsr/V_nsr);

		// Ca_JSR
		dY[22] = J_tr_j - Y[27] - J_ip3r;

		// Ca_CSR
		dY[23] = J_tr_c - Y[26];

		// Ca_PCS
		dY[24] = Beta_pcs*((2.0*I_naca_pcs - i_cal)*1000.0/(2.0*F*V_pcs*cell_volume) + (Y[27] + J_ip3r)*(V_jsr/V_pcs) - J_diff);

		// CAMK_trap
		dY[25] = alpha_camk*CAMK_bound*(CAMK_bound+Y[25])-beta_camk*Y[25];

		//Ryr2
		dY[26] = (Ryr2_inf - Y[26])/t_ryr2;

		//Ryr3
		dY[27] = (Ryr3_inf - Y[27])/t_ryr3;

		//uip3cr
		dY[28] = Y[24]*k2*(1.0-Y[28])-k2a*Y[28];

		// - Non-calcium intracellular ion concentrations -

		// K_i
		dY[3] = (-1.0*i_to_fast+-1.0*i_to_sustained+-1.0*i_kr+-1.0*i_ks+-1.0*i_k1+-1.0*i_kb+-1.0*i_f_k+-1.0*i_stim+2.0*i_nak)*1000.0/(F*diffusable_volume);

		// Na_i
		dY[4] = (-3.0*i_naca)*1000.0/(F*diffusable_volume) + (-1.0*i_na_fast+-1.0*i_na_late+-3.0*i_nak+-1.0*i_f_na+-1.0*i_nab)*1000.0/(F*diffusable_volume);

		// - x-Gates -

		// xCaL
		double x_inf_Ltype = 1.0/(1.0+exp( (Y[5]+14.6)/-5.5));
		dY[6] = (x_inf_Ltype-Y[6])/tau_x_Ltype;
		// xCaT
		double x_inf_Ttype = 1.0/(1.0+exp((Y[5]+47.8)/-5.5));
		dY[7] = (x_inf_Ttype-Y[7])/tau_x_Ttype;

		// xKs
		double x_inf_ks = 1.0/(1.0+exp((Y[5]-1.5)/-16.7));
		double tau_x_ks = 0.0;

		if (abs(Y[5]+30.0) < 0.0145)
		{
			tau_x_ks = 417.9462;
		}
		else
		{
			tau_x_ks = 1.0/(0.0000719*(Y[5]+30.0)/(1.0-exp(-0.148*(Y[5]+30.0)))+0.000131*(Y[5]+30.0)/(exp(0.0687*(Y[5]+30.0))-1.0));
		}

		dY[8] = (x_inf_ks-Y[8])/tau_x_ks;
		// xNaF
		double x_inf_na_fast = 1.0/(1.0+exp((Y[5]+25.0)/-5.0));
		dY[9] = (x_inf_na_fast-Y[9])/tau_x_na_fast;
		// xNaL
		double x_inf_na_late = 1.0/(1.0+exp((Y[5]+30.0)/-5.0));
		dY[10] = (x_inf_na_late-Y[10])/tau_x_na_late;
		// xToFast
		double x_inf_to_fast = 1.0/(1.0+exp((Y[5]+(-7.0))/-9.0));
		dY[11] = (x_inf_to_fast-Y[11])/tau_x_to_fast;

		// - y-Gates -

		// yCaL
		double y_inf_Ltype = 1.0/(1.0+exp((Y[5]+31.0)/5.54));
		double tau_y_Ltype = 25.1/(0.04+0.7*exp(-0.025*pow((Y[5]+14.5),2.0)));
		dY[12] = (y_inf_Ltype-Y[12])/tau_y_Ltype;
		// yCaT
		double y_inf_Ttype = 1.0/(1.0+exp((Y[5]+67.9)/3.87));
		double tau_y_Ttype = 1.42271*exp(-0.05119*Y[5]);
		dY[13] = (y_inf_Ttype-Y[13])/tau_y_Ttype;
		// yCaL_Ca
		double y_ca_inf_Ltype = 0.4+0.6/(1.0+pow((Y[24]/0.0001),2.0));
		double tau_y_ca_Ltype = 2.0+80.0/(1.0+pow((Y[24]/0.0001),2.0));
		dY[14] = (y_ca_inf_Ltype-Y[14])/tau_y_ca_Ltype;
		// y_f
		double y_inf_f_gate_1 = 1.0/(1.0+exp((Y[5]+109.0)/10.0));
		double tau_y_f_gate_1 = 6000.0/(exp(-1.0*(2.9+0.04*Y[5]))+exp(1.0*(3.6+0.11*Y[5])));
		dY[15] = (y_inf_f_gate_1-Y[15])/tau_y_f_gate_1;
		// y_f2
		double y_inf_f_gate_2 = 1.0/(1.0+exp((Y[5]+109.0)/10.0));
		double tau_y_f_gate_2 = 6000.0/(exp(-1.0*(2.9+0.04*Y[5]))+exp(1.0*(3.6+0.11*Y[5])));
		dY[16] = (y_inf_f_gate_2-Y[16])/tau_y_f_gate_2;

		// yKr
		double y_inf_kr = 1.0/(1.0+exp((Y[5]+50.0)/-7.5));

		double ykrv1 = 0.0;

		if (abs(Y[5]+7.0) > 0.001)
		{
			ykrv1 = 0.00138*1.0*(Y[5]+7.0)/(1.0-exp(-0.123*(Y[5]+7.0)));
		}
		else
		{
			ykrv1 = 0.00138/0.123;
		}

		double ykrv2 = 0.0;

		if (abs(Y[5]+10.0) > 0.001)
		{
			ykrv2 = 0.000061*1.0*(Y[5]+10.0)/(exp(0.145*(Y[5]+10.0))-1.0);
		}
		else
		{
			ykrv2 = 0.00061/0.145;
		}

		double tau_y_kr = 1.0/(ykrv1+ykrv2);
		dY[17] = (y_inf_kr-Y[17])/tau_y_kr;

		// yKs
		double y_inf_ks = x_inf_ks;
		double tau_y_ks = 4.0*tau_x_ks;
		dY[18] = (y_inf_ks-Y[18])/tau_y_ks;
		// yNaF
		double y_inf_na_fast = 1.0/(1.0+exp((Y[5]+69.0)/3.96));
		double tau_y_na_fast = 0.08+2.0/(1.0+exp((Y[5]+30.0)/10.0));
		// double tau_y_na_fast = 2; // Model change to improve fit of peak Vm to experiments
		dY[19] = (y_inf_na_fast-Y[19])/tau_y_na_fast;
		// yNaL
		double y_inf_na_late = 0.1+0.9/(1.0+exp((Y[5]+75.6)/6.3));
		double tau_y_na_late = 120.0+1.0*exp((Y[5]+100.0)/25.0);
		dY[20] = (y_inf_na_late-Y[20])/tau_y_na_late;
		// yToFast
		double y_inf_to_fast = 1.0/(1.0+exp((Y[5]+27.5)/8.0));
		dY[21] = (y_inf_to_fast-Y[21])/tau_y_to_fast;

		// --- Model implementation specific code ---

		//Copy dY vector into N_dY N_vector
		//CopyFromStdVector(dY, N_dY);
		for (int i=0; i<29; i++)
		{
			SetVectorComponent(N_dY,i,dY[i]);
			//SetVectorComponent(N_dY,i,-(dY[i]-i));
		}
	}
};

template<>
void OdeSystemInformation<OlliePurkinjeCvode>::Initialise()
{
	this->mVariableNames.push_back("Ca_i_bulk");
	this->mVariableUnits.push_back("Concentration");
	this->mInitialConditions.push_back(2.7012915e-05);

	this->mVariableNames.push_back("Ca_i_peripheral");
	this->mVariableUnits.push_back("Concentration");
	this->mInitialConditions.push_back(4.4031329e-05);

	this->mVariableNames.push_back("Ca_sr");
	this->mVariableUnits.push_back("Concentration");
	this->mInitialConditions.push_back(8.3420106e-01);

	this->mVariableNames.push_back("K_i");
	this->mVariableUnits.push_back("Concentration");
	this->mInitialConditions.push_back(1.1691364e+02);

	this->mVariableNames.push_back("Na_i");
	this->mVariableUnits.push_back("Concentration");
	this->mInitialConditions.push_back(5.3733057e+00);

	this->mVariableNames.push_back("Vm");
	this->mVariableUnits.push_back("mV");
	this->mInitialConditions.push_back(-8.3887820e+01);

	this->mVariableNames.push_back("x_Ltype");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(3.2082567e-06);

	this->mVariableNames.push_back("x_Ttype");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(1.3402013e-03);

	this->mVariableNames.push_back("x_ks");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(1.9477106e-02);

	this->mVariableNames.push_back("x_na_fast");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(7.5783194e-06);

	this->mVariableNames.push_back("x_na_late");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(1.9648050e-05);

	this->mVariableNames.push_back("x_to_fast");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(3.9840055e-05);

	this->mVariableNames.push_back("y_Ltype");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(7.9043210e-01);

	this->mVariableNames.push_back("y_Ttype");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(9.2061807e-01);

	this->mVariableNames.push_back("y_ca_Ltype");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(9.1681339e-01);

	this->mVariableNames.push_back("y_gate_f_k");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(4.9939260e-03);

	this->mVariableNames.push_back("y_gate_f_na");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(4.9939260e-03);

	this->mVariableNames.push_back("y_kr");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(2.3977136e-01);

	this->mVariableNames.push_back("y_ks");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(4.5066050e-02);

	this->mVariableNames.push_back("y_na_fast");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(9.7879837e-01);

	this->mVariableNames.push_back("y_na_late");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(7.4559174e-01);

	this->mVariableNames.push_back("y_to_fast");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(8.4424020e-01);

	this->mVariableNames.push_back("Ca_SSR");
	this->mVariableUnits.push_back("Concentration");
	this->mInitialConditions.push_back(8.2220728e-01);

	this->mVariableNames.push_back("Ca_CSR");
	this->mVariableUnits.push_back("Concentration");
	this->mInitialConditions.push_back(8.3724729e-01);

	this->mVariableNames.push_back("Ca_PCS");
	this->mVariableUnits.push_back("Concentration");
	this->mInitialConditions.push_back(3.8961998e-05);

	this->mVariableNames.push_back("CAMK_trap");
	this->mVariableUnits.push_back("");
	this->mInitialConditions.push_back(3.9147469e-03);

	this->mVariableNames.push_back("RyR2");
	this->mVariableUnits.push_back("");
	this->mInitialConditions.push_back(6.5158637e-25);

	this->mVariableNames.push_back("RyR3");
	this->mVariableUnits.push_back("");
	this->mInitialConditions.push_back(8.9061434e-05);

	this->mVariableNames.push_back("UIP3CR");
	this->mVariableUnits.push_back("");
	this->mInitialConditions.push_back(2.5035068e-01);

	this->mInitialised = true;
}

#endif
