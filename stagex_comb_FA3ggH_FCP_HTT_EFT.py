from HiggsAnalysis.CombinedLimit.PhysicsModel import *
import fnmatch 

class stagex_comb(PhysicsModel):
    "Allow different signal strength fits for the stage-x model"

    xsecs = {
            "sigma1_HZZ": 290.58626,
            #          "sigma3_HZZ": 581.17253,
            "sigma3_HZZ": 44.670158,
            "sigma1_VBF": 968.674,
            "sigma3_VBF": 10909.54,
            "sigma1_ZH":  9022.36,
            "sigma3_ZH":  434763.7,
            "sigma1_WH":  30998.54,
            "sigma3_WH":  2028656,
            "sigma2_HZZ": 105.85594,
            "sigmaL1_HZZ": 1.9846071e-06,
            "sigma2_VBF": 13102.71,
            "sigmaL1_VBF": 2.08309E-4,
            "sigma2_ZH": 713123,
            "sigmaL1_ZH": 33652.46e-6,
            "sigma2_WH": 3106339,
            "sigmaL1_WH": 11234.91e-5,
            "sigmaa1a3int_VBF": 1937.15,
            "sigmaa1a3int_ZH": 18044.72,
            "sigmaa1a3int_WH": 61997.07,
            # for ggH what is listed are yields in histos (all normliazed to powheg xsection):
            "yield_Powheg_ggH": 5.127128e+00,
            "yield_SM_ggH": 7.937751e+01,
            "yield_BSM_ggH": 8.952848e+01,
            "sigma_Powheg_ggH": 48.58,
            "sigma_SM_ggH": 15980,
            "sigma_BSM_ggH": 15981,
            "BR_H_tt": 0.0627

        }

    def __init__(self):
        self.POIs = ""
        self.actth=False 
        self.acggh=False 
        self.accor=False 
        self.useint=False 
        self.eft=False 
        self.options= ""
        self.stage0= False 
        self.rvrf= False 
        self.singlemu= False 
        self.splitHad= False 
        self.hggalone= False 
        self.muNames =[]
        self.pois = []
        self.count=0
        self.doHTTHZZ= False 
    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self, modelBuilder)
        self.modelBuilder.doModelBOnly = False
    def getYieldScale(self,bin,process):
            if process in ["GGH2Jets_sm_M",]:
                return 'smCoupling_ggH'
            if process in ["GGH2Jets_pseudoscalar_M",]:
                return 'bsmCoupling_ggH'
            if process in ["reweighted_qqH_htt_0PM",]:
                return 'smCoupling_VBF'
            if process in ["reweighted_WH_htt_0PM",]:
                return 'smCoupling_WH'
            if process in ["reweighted_ZH_htt_0PM",]:
                return 'smCoupling_ZH'
            if process in ["reweighted_qqH_htt_0M",]:
                return 'bsmCoupling_VBF'
            if process in ["reweighted_WH_htt_0M",]:
                return 'bsmCoupling_WH'
            if process in ["reweighted_ZH_htt_0M",]:
                return 'bsmCoupling_ZH'
            if process in ["reweighted_qqH_htt_0Mf05ph0"]:
                return 'intCoupling_VBF'
            if process in ["reweighted_ZH_htt_0Mf05ph0"]:
                return 'intCoupling_ZH'
            if process in ["reweighted_WH_htt_0Mf05ph0"]:
                return 'intCoupling_WH'
            if process in ["GGH2Jets_pseudoscalar_Mf05ph0"]:
                return 'intCoupling_ggH'

            print "",process, self.DC.isSignal[process]
            if not self.DC.isSignal[process]: 
                    return 1
            else:
                muname=""
                if self.stage0 :
                        if process.startswith("ggH") or process.startswith("ggh"):
                                if self.acggh: 
                                    if "ALT" in process:
                                        if self.eft:
                                                muname= "c_ts"
                                        else:
                                                muname= "r_ggH_times_x"
                                    elif "Mix" in process:
                                        if self.eft:
                                                if self.useint:
                                                        muname= "c_int"
                                                else:
                                                        return 0
                                        else:
                                                if self.useint:
                                                        muname= "r_ggH_times_int"
                                                else:
                                                        return 0
                                    else:
                                        if self.eft:
                                                if "125" in process:
                                                        muname = "c_s_all"
                                                else:
                                                        muname= "c_s"
                                        else:
                                                muname= "r_ggH_times_notx"
                                else:
                                    if not "ALT" in process and not "Mix" in process:
                                        if self.eft:
                                                return 1
                                        else:
                                           muname = "r_ggH"
                                    else:
                                            return 0
                        elif "VHori" in process :
                                muname= "r_VH"
                        elif "VBFori" in process :
                                muname= "r_VBF"
                        elif process.startswith("VBF") :
                                if ("ALT" in process):
                                  muname= "gsVBF_ALT"
                                elif ("Mix" in process):
                                  muname= "gsVBF_mix"
                                else:
                                  muname = "gsVBF"
                        elif process.startswith("WH"): 
                                muname= "gsWH"
                        elif process.startswith("ZH"):
                                muname= "gsZH"
                        elif process.startswith("tqH") or "thq" in process or "thw" in process:
                             if self.actth:
                                if self.eft:
                                   muname="kappa_th"
                                else: 
                                   muname="r_TH_times_ac"
                             else:
                                muname="r_TTH"
                        elif process.startswith("TTH") or process.startswith("tth"): 
                                if self.actth:
                                    if "ALT" in process:
                                        if self.eft:
                                                muname= "kappa_ts"
                                        else:
                                                muname= "r_TTH_times_x"
                                    else:
                                        if self.eft:
                                                muname= "kappa_s"
                                        else:
                                                muname= "r_TTH_times_notx"
                                else:
                                    if not "ALT" in process:
                                        if self.eft:
                                                return 1
                                        else:
                                           muname = "r_TTH"
                                    else:
                                        return 0
#                       elif process.startswith("TH") or process.startswith("tqH"):
#                               muname= "r_TTH"
                        elif process.startswith("BBH"):
                                muname= "r_ggH"
                        elif process.startswith("VH"):
                                muname= "g1s"
                        if self.splitHad:       
                                if ("VHori" in process):
                                        muname= "r_VH_Had"
                                elif (process.startswith("VH") and (not "VBFori" in process)):
                                        muname= "r_VH_Had"
                                if "Lep" in process:
                                        muname= "r_VH_Lep"
                        if "125" in process and not self.hggalone:
                                muname += "_times_rhgg"
                elif self.singlemu:
                        muname = "r"
                elif self.rvrf:
                        if ("VH" in process or "VBF" in process):
                                muname= "rv"
                        elif ("ggH" in process or "BBH" in process or "TH" in process):
                                muname= "rf"
                else:
                        muname = "r_%s"%process
                        if process.startswith("BBH"):
                                muname= "r_ggH_0j_10_200"
                        if process.startswith("TH"):
                                muname= "r_TTH"
                        muname = muname.replace("_VHori","") 
                        muname = muname.replace("_VBFori","") 
                if self.modelBuilder.out.var(muname):
                        print "reclying %s" %muname
                else:
                        if muname.startswith("r_") and not "times" in muname:
                                self.modelBuilder.doVar("%s[1,0,10]" % muname)
                                print "scale process %s with %s" %(process,muname)
                                self.pois.append(muname)
                                self.POIs=",".join(self.pois)
                                self.modelBuilder.doSet("POI",self.POIs)
                                print "Default parameters of interest: ", self.POIs
                print "final muname ", muname
                return muname 
    def setPhysicsOptions(self,physOptions):
            for po in physOptions:
                    if 'doStage0' in po: 
                            self.stage0= True
                            print "doing stage0"
                    if 'doactth' in po: 
                            self.actth= True
                            print "doing tth AC ttH"
                    if 'doaccor' in po: 
                            self.accor= True
                            print "doing tth AC ttH"
                    if 'doacggh' in po: 
                            self.acggh= True
                            print "doing tth AC ggH"
                    if 'useint' in po: 
                            self.useint= True
                            print "use interference"
                    if 'hggalone' in po: 
                            self.hggalone= True
                            print "hggalone"
                    if 'doeft' in po: 
                            self.eft= True
                            print "doing EFT framework"
                    if 'singlemu' in po: 
                            self.singlemu= True
                            print "doing single mu"
                    if 'rvrf' in po: 
                            self.rvrf= True
                            print "doing rvrf"
                    if 'splitHad' in po: 
                            self.splitHad= True
                            print "Splitting had and lep VH"
                    if 'doHTTHZZ' in po: 
                            self.doHTTHZZ= True
                            print "doing HTT+HZZ"

    def doParametersOfInterest(self):
            self.POIs=",".join(self.pois)
            if self.actth or self.acggh:
                #self.modelBuilder.doVar("g1s[1,0,10]" )
                #self.modelBuilder.factory_("expr::gHz(\"(@0*@0)\", a1)");
#               self.modelBuilder.doVar("gHz[1]");

###This is the original one
                self.modelBuilder.doVar("a1[1,0,10]" )
                self.modelBuilder.doVar("a3[0,-1,1]" )
                self.modelBuilder.factory_("expr::gHz(\"(@0*@0+@1*@1*0.1537)\", a1,a3)");
                ########
                #######Simpler version
                #self.modelBuilder.doVar("rV[1,0,10]" )
                #self.modelBuilder.doVar("a3[0,-1,1]" )
                #self.modelBuilder.factory_("expr::gsZH(\"@0\", rV)");
                #self.modelBuilder.factory_("expr::gsWH(\"@0\", rV)");

                #self.modelBuilder.factory_("expr::gsVBF(\"@0*(1-@1*@1)\", rV,a3)");
                #self.modelBuilder.factory_("expr::gsVBF_ALT(\"@0*(@1*@1)\", rV,a3)");
                #self.modelBuilder.factory_("expr::gsVBF_mix(\"0\", rV)");

                #self.modelBuilder.factory_("expr::gsVBF(\"(@0*@0-3.356*@0*@1)*@2\", a1,a3,gHz)");
                #self.modelBuilder.factory_("expr::gsVBF_ALT(\"(@0*@0*11.2627-3.356*@0*@1)*@2\", a3,a1,gHz)");
                #self.modelBuilder.factory_("expr::gsVBF_mix(\"(@0*@1*3.356*2*@2)\", a1,a3,gHz)");
                #self.modelBuilder.factory_("expr::gsVBF(\"(@0*@0)*@2\", a1,a3,gHz)");
                #self.modelBuilder.factory_("expr::gsVBF_ALT(\"(@0*@0*11.2627)*@2\", a3,a1,gHz)");
                #self.modelBuilder.factory_("expr::gsVBF_mix(\"0\", a1,a3,gHz)");
                if self.eft:
                        self.modelBuilder.doVar("kappa[1,0,10]" )
                        self.modelBuilder.doVar("kappa_t[0,-10,10]" )
                        self.modelBuilder.doVar("cgg[0,-10,10]" )
                        self.modelBuilder.doVar("cgg_t[0,-10,10]" )
                        
                        self.modelBuilder.doVar("kappa_tau[1,0,10]" )
                        if self.doHTTHZZ:
                                self.modelBuilder.out.var("kappa_tau").setConstant(False)
                                self.modelBuilder.out.var("kappa_tau").setAttribute("flatParam")
                        else:
                                print("Fixing kappa_tau to 1")
                                self.modelBuilder.out.var("kappa_tau").setConstant(True)
                        #125
                        #self.modelBuilder.factory_("expr::gamma(\"0.5824*(@0*@0+@1*@1)+(@3*@3*0.323+@2*@2)*0.2137+0.02619*(@2*@2+0.152*@3*@3)+0.08187*(@1*@1*2.383+@0*@0)+0.06272+0.02891+0.0002176+0.00227*(1.6057+0.07135*@0*@0-0.6770*@0*sqrt(@2)+0.1668*@1*@1)+0.001533*(1.1183+0.0033*@0*@0+0.01213*@1*@1-0.1219*@0*sqrt(@2))\", kappa,kappa_t,a1,a3)");
                        #Hbb,HWW+HZZ,Hgg,Htautau,Hcc,Hmm,Hgammagamma,Hzg 
                        #125.38
                        #self.modelBuilder.factory_("expr::gamma(\"0.5760*(@0*@0+@1*@1)+(@3*@3*0.323+@2*@2)*0.2203+0.02716*(@2*@2+0.152*@3*@3)+0.08154*(@1*@1*2.382+@0*@0+@4*@4*1.0298-1.3204*@0*@4+1.0298*@5*@5+3.1294*@5*@1)+0.06208+0.02860+0.0002153+0.00227*(1.6057+0.07135*@0*@0-0.6770*@0*sqrt(@2)+0.1668*@1*@1)+0.001567*(1.1183+0.0033*@0*@0+0.01213*@1*@1-0.1219*@0*sqrt(@2))\", kappa,kappa_t,a1,a3,cgg,cgg_t)");
                        if not self.accor:
                                if self.acggh:
#                                               self.modelBuilder.factory_("expr::gamma(\"0.5824+(@3*@3*0.323+@2*@2)*0.2137+0.02619*(@2*@2+0.152*@3*@3)+0.08187*(@1*@1+@0*@0)+0.06272+0.02891+0.0002176+0.00227+0.001533\", cgg,cgg_t,a1,a3)"); 
                                        # self.modelBuilder.factory_("expr::gamma(\"0.5760+(@3*@3*0.323+@2*@2)*0.2203+0.02716*(@2*@2+0.152*@3*@3)+0.08154*(1+@0*@0*1.0298+2.0297*@0+1.0298*@1*@1)+0.06208+0.02860+0.0002153+0.00227+0.001567\", cgg,cgg_t,a1,a3)");
                                        self.modelBuilder.factory_("expr::gamma(\"0.5760+(@3*@3*0.323+@2*@2)*0.2203+0.02716*(@2*@2+0.152*@3*@3)+0.08154*(1+@0*@0*1.0298+2.0297*@0+1.0298*@1*@1)+0.06208*@4*@4+0.02860+0.0002153+0.00227+0.001567\", cgg,cgg_t,a1,a3,kappa_tau)");
                                else:
                                        self.modelBuilder.factory_("expr::gamma(\"0.5760+(@3*@3*0.323+@2*@2)*0.2203+0.02716*(@2*@2+0.152*@3*@3)+0.08154+0.06208+0.02860+0.0002153+0.00227+0.001567\", kappa,kappa_t,a1,a3)");
                                if self.useint:
                                        self.modelBuilder.factory_("expr::c_ts(\"(@0*@0*1.0298+1.0298*@0*@1+1.0148*@0)/@2\", cgg_t,cgg,gamma)");
                                        self.modelBuilder.factory_("expr::c_s(\"(1+@0*@0*1.0298+2.0297*@0+1.0298*@0*@1+1.0148*@1)/@2\", cgg,cgg_t,gamma)");
                                else:
                                        self.modelBuilder.factory_("expr::c_ts(\"@0*@0*1.0298/@1\", cgg_t,gamma)");
                                        self.modelBuilder.factory_("expr::c_s(\"(1+@0*@0*1.0298+2.0297*@0)/@1\", cgg,gamma)");
                                self.modelBuilder.factory_("expr::c_int(\"-2*(1.0298*@0*@1+1.0148*@1)/@2\", cgg, cgg_t,gamma)");
                                self.modelBuilder.factory_("expr::kappa_hgg(\"1\", kappa)");
                        else:
#                                       self.modelBuilder.factory_("expr::gamma(\"0.5760*(@0*@0+@1*@1)+(@3*@3*0.323+@2*@2)*0.2203+0.02716*(@2*@2+0.152*@3*@3)+0.08154*(@1*@1*2.382+@0*@0+@4*@4*1.0298-1.3204*@0*@4+1.0298*@5*@5+3.1294*@5*@1)+0.06208+0.02860+0.0002153+0.00227*(1.6057+0.07135*@0*@0-0.6770*@0*sqrt(@2)+0.1668*@1*@1)+0.001567*(1.1183+0.0033*@0*@0+0.01213*@1*@1-0.1219*@0*sqrt(@2))\", kappa,kappa_t,a1,a3,cgg,cgg_t)");
                                # self.modelBuilder.factory_("expr::gamma(\"0.5760+(@3*@3*0.323+@2*@2)*0.2203+0.02716*(@2*@2+0.152*@3*@3)+0.08154*(pow((1.0520 *@0 + 1.0148 *@4),2) + pow((1.6036*@1+1.0148*@5),2) - 0.1093 *(1.0520 *@0 + 1.0148 *@4)  + 0.0082)+0.06208+0.02860+0.0002153+0.00227*(1.6057*@2*@2+0.07132*@0*@0-0.6855*@0*@2+0.1699*@1*@1-0.0018*@0+0.0085*@2+0.00002)+0.001567*(1.1186*@2*@2+0.0035*@0*@0+0.0126*@1*@1-0.125*@0*@2+0.000003-0.000183*@0+ 0.0031*@2)\", kappa,kappa_t,a1,a3,cgg,cgg_t)");
                                self.modelBuilder.factory_("expr::gamma(\"0.5760+(@3*@3*0.323+@2*@2)*0.2203+0.02716*(@2*@2+0.152*@3*@3)+0.08154*(pow((1.0520 *@0 + 1.0148 *@4),2) + pow((1.6036*@1+1.0148*@5),2) - 0.1093 *(1.0520 *@0 + 1.0148 *@4)  + 0.0082)+0.06208*@6*@6+0.02860+0.0002153+0.00227*(1.6057*@2*@2+0.07132*@0*@0-0.6855*@0*@2+0.1699*@1*@1-0.0018*@0+0.0085*@2+0.00002)+0.001567*(1.1186*@2*@2+0.0035*@0*@0+0.0126*@1*@1-0.125*@0*@2+0.000003-0.000183*@0+ 0.0031*@2)\", kappa,kappa_t,a1,a3,cgg,cgg_t,kappa_tau)");
#                               ( 1.0520 t + 1.0148 c)^2 + (1.6036t’+1.0148c’)^2 - 0.1093 *(1.0520 t + 1.0148 c)  + 0.0082
#  pow((1.0520 *@0 + 1.0148 *@4),2) + pow((1.6036*@1+1.0148*@5),2) - 0.1093 *(1.0520 *@0 + 1.0148 *@4)  + 0.0082
#                               @1*@1*2.5717+1.1068*@0*@0+@4*@4*1.0298+2.1357*@0*@4+1.0298*@5*@5+3.2547*@5*@1-0.1150*@0-0.1109*@4+0.0082

                                if self.useint:
                                        self.modelBuilder.factory_("expr::c_ts(\"(pow((1.6036*@0 + 1.0148*@1),2)+ (1.6036*@0 + 1.0148*@1)*(1.0520*@2 + 1.0148*@3))/@4\", kappa_t,cgg_t,kappa,cgg,gamma)");
                                        self.modelBuilder.factory_("expr::c_s(\"( pow((1.0520*@0 + 1.0148*@1),2)- 0.1093 *(1.0520 *@0 + 1.0148 *@1) +0.0082+(1.6036*@2 + 1.0148*@3)*(1.0520*@0 + 1.0148*@1))/@4\", kappa,cgg,kappa_t,cgg_t,gamma)");
                                else:
                                        self.modelBuilder.factory_("expr::c_ts(\"(pow((1.6036*@0 + 1.0148*@1),2))/@2\", kappa_t,cgg_t,gamma)");
                                        self.modelBuilder.factory_("expr::c_s(\"(pow((1.0520*@0 + 1.0148*@1),2)- 0.1093 *(1.0520 *@0 + 1.0148 *@1) +0.0082)/@2\", kappa,cgg,gamma)");
#                                       self.modelBuilder.factory_("expr::kappa_hgg(\"1.604+0.0681*@0*@0-0.673*@0+0.1611*@1*@1\", kappa,kappa_t)");
                                self.modelBuilder.factory_("expr::kappa_hgg(\"1.6057*@2*@2+0.07132*@0*@0-0.6855*@0*@2+0.1699*@1*@1-0.0018*@0+0.0085*@2+0.00002\", kappa,kappa_t,a1)");
                                self.modelBuilder.factory_("expr::c_s_all_times_rhgg(\"(@0*@0+@1*@1*2.382)*@2/@3\", kappa,kappa_t,kappa_hgg,gamma)");
                                self.modelBuilder.factory_("expr::c_int(\"(-2*(1.6036*@0 + 1.0148*@1)*(1.0520*@2 + 1.0148*@3))/@4\",kappa_t, cgg_t,kappa,cgg, gamma)");
                                self.modelBuilder.factory_("expr::c_ts_times_rhgg(\"@0*@1\", c_ts,kappa_hgg)");
                                self.modelBuilder.factory_("expr::c_s_times_rhgg(\"@0*@1\", c_s, kappa_hgg)");

                        self.modelBuilder.factory_("expr::kappa_ts(\"@0*@0/2.56/@1\", kappa_t,gamma)");
                        self.modelBuilder.factory_("expr::kappa_s(\"@0*@0/@1\", kappa,gamma)");
                        self.modelBuilder.factory_("expr::kappa_th(\"(3.067*@0*@0+3.590*@2*@2-5.657*@0*@2+0.975*@1*@1)/@3\", kappa,kappa_t,a1,gamma)");
       
                        self.modelBuilder.factory_("expr::kappa_ts_times_rhgg(\"@0*@1\",kappa_hgg,kappa_ts)");
                        self.modelBuilder.factory_("expr::kappa_s_times_rhgg(\"@0*@1\", kappa_hgg,kappa_s)");
                        self.modelBuilder.factory_("expr::kappa_th_times_rhgg(\"@0*@1\", kappa_hgg,kappa_th)");

                        self.modelBuilder.factory_("expr::gsZH(\"(@0*@0+@1*@1*48.19)*@2/@3\", a1,a3,gHz,gamma)");
                        self.modelBuilder.factory_("expr::gsWH(\"(@0*@0+@1*@1*65.44)*@2/@3\", a1,a3,gHz,gamma)");

                        self.modelBuilder.factory_("expr::gsVBF(\"(@0*@0-@0*@1)*@2/@3\", a1,a3,gHz,gamma)");
                        self.modelBuilder.factory_("expr::gsVBF_ALT(\"(@0*@0-@0*@1)*@2/@3\", a3,a1,gHz,gamma)");
                        self.modelBuilder.factory_("expr::gsVBF_mix(\"(@0*@1*@2)/@3\", a1,a3,gHz,gamma)");
       
                        self.modelBuilder.factory_("expr::bsmCoupling_ggH(\"@0*(@1**2)\", c_ts, kappa_tau)");
                        self.modelBuilder.factory_("expr::smCoupling_ggH(\"@0*(@1**2)\", c_s, kappa_tau)");
                        self.modelBuilder.factory_("expr::intCoupling_ggH(\"@0*(@1**2)\", c_int, kappa_tau)");

                        self.modelBuilder.doVar('expr::muV("(@0**2)*(@1**2 + @2**2*9.772)", kappa_tau, a1, a3)')
                        self.modelBuilder.doVar('expr::CMS_zz4l_fai1("@1**2*0.1537/(@1**2*0.1537+@0**2)", a1, a3)')
                        self.modelBuilder.doVar('expr::a1_HTT("sqrt(1-abs(@0))", CMS_zz4l_fai1)')
                        self.modelBuilder.doVar('expr::a3_HTT("(@0>0 ? 1 : -1) * sqrt(abs(@0)*{sigma1_HZZ}/{sigma3_HZZ})", CMS_zz4l_fai1)'.format(**stagex_comb.xsecs))

                        self.modelBuilder.factory_('expr::smCoupling_VBF("@0*@1**2 - @0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})", muV,a1_HTT,a3_HTT)'.format(**stagex_comb.xsecs))
                        self.modelBuilder.factory_('expr::bsmCoupling_VBF("@0*@1**2*{sigma3_VBF}/{sigma1_VBF} - @0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})", muV,a3_HTT,a1_HTT)'.format(**stagex_comb.xsecs))
                        self.modelBuilder.factory_('expr::intCoupling_VBF("@0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})*{sigmaa1a3int_VBF}/{sigma1_VBF}", muV,a1_HTT,a3_HTT)'.format(**stagex_comb.xsecs))
                        
                        self.modelBuilder.factory_('expr::smCoupling_ZH("@0*@1**2 - @0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})", muV,a1_HTT,a3_HTT)'.format(**stagex_comb.xsecs))
                        self.modelBuilder.factory_('expr::bsmCoupling_ZH("@0*@1**2*{sigma3_ZH}/{sigma1_ZH} - @0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})", muV,a3_HTT,a1_HTT)'.format(**stagex_comb.xsecs))
                        self.modelBuilder.factory_('expr::intCoupling_ZH("@0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})*{sigmaa1a3int_ZH}/{sigma1_ZH}", muV,a1_HTT,a3_HTT)'.format(**stagex_comb.xsecs))
                        
                        self.modelBuilder.factory_('expr::smCoupling_WH("@0*@1**2 - @0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})", muV,a1_HTT,a3_HTT)'.format(**stagex_comb.xsecs))
                        self.modelBuilder.factory_('expr::bsmCoupling_WH("@0*@1**2*{sigma3_WH}/{sigma1_WH} - @0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})", muV,a3_HTT,a1_HTT)'.format(**stagex_comb.xsecs))
                        self.modelBuilder.factory_('expr::intCoupling_WH("@0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})*{sigmaa1a3int_WH}/{sigma1_WH}", muV,a1_HTT,a3_HTT)'.format(**stagex_comb.xsecs))

                        if self.actth:
                                self.pois.append("kappa,kappa_t")
                        if self.acggh:
                                self.pois.append("cgg,cgg_t")
       
                else:
                        self.modelBuilder.doVar("r_TTH[1,0,10]" )
                        self.modelBuilder.doVar("r_ggH[1,0,10]" )
#                               self.modelBuilder.doVar("rhgg[1,0,10]" )
                        self.modelBuilder.doVar("r_TTH_gg[1,0,10]" )
                        self.modelBuilder.doVar("x[0,-1,1]" )
                        self.modelBuilder.factory_("expr::xinmu(\"@0/2.56/(1-abs(@0)+abs(@0)/2.56)\", x)");
#                               self.modelBuilder.factory_("expr::xinmuggh(\"@0*2.381/(1-abs(@0)+abs(@0)*2.381)\", x)");
                        self.modelBuilder.factory_("expr::r_TTH_times_x(\"@0*abs(@1)\", r_TTH, xinmu)");
                        self.modelBuilder.factory_("expr::r_TTH_times_notx(\"@0*(1-abs(@1))\", r_TTH, xinmu)");
                        self.modelBuilder.factory_("expr::r_TH_times_ac(\"3.067*@0*(1-abs(@1))+3.590-5.657*sqrt(@0*(1-abs(@1)))+0.975*abs(@1)*@0*2.56\",  r_TTH,xinmu)");
                        self.modelBuilder.factory_("expr::r_TTH_times_x_times_rhgg(\"@0*abs(@1)\", r_TTH_gg,xinmu)");
                        self.modelBuilder.factory_("expr::r_TTH_times_notx_times_rhgg(\"@0*(1-abs(@1))\", r_TTH_gg,xinmu)");
                        self.modelBuilder.factory_("expr::r_TH_times_ac_times_rhgg(\"3.067*@0*(1-abs(@1))+3.590-5.657*sqrt(@0*(1-abs(@1)))+0.975*abs(@1)*@0*2.56\",  r_TTH_gg,xinmu)");
                        self.modelBuilder.factory_("expr::gsZH(\"(@0*@0+@1*@1*48.19)*@2\", a1,a3,gHz)");
                        self.modelBuilder.factory_("expr::gsWH(\"(@0*@0+@1*@1*65.44)*@2\", a1,a3,gHz)");

                        self.modelBuilder.factory_("expr::gsVBF(\"(@0*@0-@0*@1)*@2\", a1,a3,gHz)");
                        self.modelBuilder.factory_("expr::gsVBF_ALT(\"(@0*@0-@0*@1)*@2\", a3,a1,gHz)");
                        self.modelBuilder.factory_("expr::gsVBF_mix(\"(@0*@1*@2)\", a1,a3,gHz)");


                        self.modelBuilder.doVar("HTT_over_HZZ[1.0,0,10]")
                        if self.doHTTHZZ:
                                self.modelBuilder.out.var("HTT_over_HZZ").setConstant(False)
                                self.modelBuilder.out.var("HTT_over_HZZ").setAttribute("flatParam")
                        else:
                                print("Fixing HTT_over_HZZ to 1")
                                self.modelBuilder.out.var("HTT_over_HZZ").setConstant(True)


                        self.modelBuilder.doVar('expr::muV("@0*@1**2 + @0*@2**2*9.772", HTT_over_HZZ, a1, a3)')
                        self.modelBuilder.doVar('expr::CMS_zz4l_fai1("@1**2*0.1537/(@1**2*0.1537+@0**2)", a1, a3)')
                        
                        self.modelBuilder.doVar('expr::a1_HTT("sqrt(1-abs(@0))", CMS_zz4l_fai1)')
                        self.modelBuilder.doVar('expr::a3_HTT("(@0>0 ? 1 : -1) * sqrt(abs(@0)*{sigma1_HZZ}/{sigma3_HZZ})", CMS_zz4l_fai1)'.format(**stagex_comb.xsecs))

                        if self.accor:
                            self.modelBuilder.doVar('expr::r_TTH_TT("@0*@1", r_TTH, HTT_over_HZZ)')
                                        
                            self.modelBuilder.factory_("expr::smCoupling_ggH(\"@0*(1-abs(@1))-1.6*(@1>0?-1.543:1.543)*@0*sqrt(1-abs(@1))*sqrt(abs(@1))\", r_TTH_TT, xinmu)")
                            self.modelBuilder.factory_("expr::bsmCoupling_ggH(\"2.382*@0*abs(@1)*2.56-1.6*(@1>0?-1.543:1.543)*@0*sqrt(1-abs(@1))*sqrt(abs(@1))\", r_TTH_TT, xinmu)")
                            self.modelBuilder.factory_("expr::intCoupling_ggH(\"1.6*(@1>0?-1.543:1.543)*2*@0*sqrt(1-abs(@1))*sqrt(abs(@1))\", r_TTH_TT, xinmu)")
                        else:
                            self.modelBuilder.doVar('expr::muf("@0*@1", r_ggH, HTT_over_HZZ)')
                            self.modelBuilder.doVar('expr::a1_ggH("sqrt(1-abs(@0))", x)')
                            self.modelBuilder.doVar('expr::a3_ggH("(@0>0 ? -1 : 1) * sqrt(abs(@0))", x)')

                            self.modelBuilder.factory_('expr::smCoupling_ggH("@0*@1**2 -  @0*@1*@2*sqrt({sigma_BSM_ggH}/{sigma_SM_ggH})", muf,a1_ggH,a3_ggH)'.format(**stagex_comb.xsecs))
                            self.modelBuilder.factory_('expr::bsmCoupling_ggH("@0*@1**2*{sigma_BSM_ggH}/{sigma_SM_ggH} - @0*@1*@2*sqrt({sigma_BSM_ggH}/{sigma_SM_ggH})", muf,a3_ggH,a1_ggH)'.format(**stagex_comb.xsecs))
                            self.modelBuilder.factory_('expr::intCoupling_ggH("@0*@1*@2*sqrt({sigma_BSM_ggH}/{sigma_SM_ggH})*2.", muf,a1_ggH,a3_ggH)'.format(**stagex_comb.xsecs))

                        self.modelBuilder.factory_('expr::smCoupling_VBF("@0*@1**2 - @0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})", muV,a1_HTT,a3_HTT)'.format(**stagex_comb.xsecs))
                        self.modelBuilder.factory_('expr::bsmCoupling_VBF("@0*@1**2*{sigma3_VBF}/{sigma1_VBF} - @0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})", muV,a3_HTT,a1_HTT)'.format(**stagex_comb.xsecs))
                        self.modelBuilder.factory_('expr::intCoupling_VBF("@0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})*{sigmaa1a3int_VBF}/{sigma1_VBF}", muV,a1_HTT,a3_HTT)'.format(**stagex_comb.xsecs))
                        
                        self.modelBuilder.factory_('expr::smCoupling_ZH("@0*@1**2 - @0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})", muV,a1_HTT,a3_HTT)'.format(**stagex_comb.xsecs))
                        self.modelBuilder.factory_('expr::bsmCoupling_ZH("@0*@1**2*{sigma3_ZH}/{sigma1_ZH} - @0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})", muV,a3_HTT,a1_HTT)'.format(**stagex_comb.xsecs))
                        self.modelBuilder.factory_('expr::intCoupling_ZH("@0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})*{sigmaa1a3int_ZH}/{sigma1_ZH}", muV,a1_HTT,a3_HTT)'.format(**stagex_comb.xsecs))
                        
                        self.modelBuilder.factory_('expr::smCoupling_WH("@0*@1**2 - @0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})", muV,a1_HTT,a3_HTT)'.format(**stagex_comb.xsecs))
                        self.modelBuilder.factory_('expr::bsmCoupling_WH("@0*@1**2*{sigma3_WH}/{sigma1_WH} - @0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})", muV,a3_HTT,a1_HTT)'.format(**stagex_comb.xsecs))
                        self.modelBuilder.factory_('expr::intCoupling_WH("@0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})*{sigmaa1a3int_WH}/{sigma1_WH}", muV,a1_HTT,a3_HTT)'.format(**stagex_comb.xsecs))

                #       self.modelBuilder.factory_("expr::r_TTH_times_x_times_rhgg(\"@0*abs(@1)*@2\", r_TTH,xinmu, rhgg)");
                #       self.modelBuilder.factory_("expr::r_TTH_times_notx_times_rhgg(\"@0*(1-abs(@1))*@2\", r_TTH,xinmu, rhgg)");
                #       self.modelBuilder.factory_("expr::r_TH_times_ac_times_rhgg(\"(3.067*@0*(1-abs(@1))+3.590-5.657*sqrt(@0*(1-abs(@1)))+0.975*abs(@1)*@0*2.56)*@2\",  r_TTH,xinmu,rhgg)");
                        if self.accor:
                                if self.useint:
                                        self.modelBuilder.factory_("expr::r_ggH_times_x(\"2.382*@0*abs(@1)*2.56-1.6*(@1>0?-1.543:1.543)*@0*sqrt(1-abs(@1))*sqrt(abs(@1))\", r_TTH, xinmu)");
                                        self.modelBuilder.factory_("expr::r_ggH_times_notx(\"@0*(1-abs(@1))-1.6*(@1>0?-1.543:1.543)*@0*sqrt(1-abs(@1))*sqrt(abs(@1))\", r_TTH, xinmu)");
                                else:
                                        self.modelBuilder.factory_("expr::r_ggH_times_x(\"2.382*@0*abs(@1)*2.56\", r_TTH, xinmu)");
                                        self.modelBuilder.factory_("expr::r_ggH_times_notx(\"@0*(1-abs(@1))\", r_TTH, xinmu)");
                                self.modelBuilder.factory_("expr::r_ggH_times_int(\"1.6*(@1>0?-1.543:1.543)*2*@0*sqrt(1-abs(@1))*sqrt(abs(@1))\", r_TTH, xinmu)");
                                self.pois.append("x,r_TTH")
                        else:
                                if self.useint:
                                        self.modelBuilder.factory_("expr::r_ggH_times_x(\"@0*abs(@1)-(@1>0?-1.0:1.0)*@0*sqrt(1-abs(@1))*sqrt(abs(@1))\", r_ggH, x)");
                                        self.modelBuilder.factory_("expr::r_ggH_times_notx(\"@0*(1-abs(@1))-(@1>0?-1.0:1.0)*@0*sqrt(1-abs(@1))*sqrt(abs(@1))\", r_ggH, x)");
                                else:
                                        self.modelBuilder.factory_("expr::r_ggH_times_x(\"@0*abs(@1)\", r_ggH, x)");
                                        self.modelBuilder.factory_("expr::r_ggH_times_notx(\"@0*(1-abs(@1))\", r_ggH, x)");
                                self.modelBuilder.factory_("expr::r_ggH_times_int(\"(@1>0?-1.0:1.0)*2*@0*sqrt(1-abs(@1))*sqrt(abs(@1))\", r_ggH, x)");
                                self.pois.append("x,r_ggH,r_TTH,r_TTH_gg")
            self.POIs=",".join(self.pois)
            print "Default parameters of interest: ", self.POIs
            self.modelBuilder.doSet("POI",self.POIs)

stagex_comb = stagex_comb()
