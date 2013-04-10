#!/bin/tcsh


echo 'Running signal samples....'

#root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(9110,"2012_53x",2,-1,1,1,false,true,false)' > & ! ttH_110_lep.log 

#echo 'DONE M=110'

#root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(9115,"2012_53x",2,-1,1,1,false,true,false)' > & ! ttH_115_lep.log 

#echo 'DONE M=115'

#root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(9120,"2012_53x",2,-1,1,1,false,true,false)' > & ! ttH_120_lep.log

#echo 'DONE M=120'

#root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(9122,"2012_53x",2,-1,1,1,false,true,false)' > & ! ttH_122_lep.log 

#echo 'DONE M=122'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(9125,"2012_53x",2,-1,1,1,false,true,false)' > & ! ttH_125_lep.log 

echo 'DONE M=125'

#root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(9127,"2012_53x",2,-1,1,1,false,true,false)' > & ! ttH_127_lep.log 

#echo 'DONE M=127'

#root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(9130,"2012_53x",2,-1,1,1,false,true,false)' > & ! ttH_130_lep.log 

#echo 'DONE M=130'

#root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(9135,"2012_53x",2,-1,1,1,false,true,false)' > & ! ttH_135_lep.log 

#echo 'DONE M=135'

#root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(9140,"2012_53x",2,-1,1,1,false,true,false)' > & ! ttH_140_lep.log 

#echo 'DONE M=140'



echo 'Running single top samples....'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2600,"2012_53x",2,-1,1,1,false,true,false)' > & ! T_s-channel_lep.log 

echo 'DONE single top s-channel'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2601,"2012_53x",2,-1,1,1,false,true,false)' > & ! Tbar_s-channel_lep.log 

echo 'DONE single top-bar s-channel'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2602,"2012_53x",2,-1,1,1,false,true,false)' > & ! T_t-channel_lep.log 

echo 'DONE single top t-channel'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2603,"2012_53x",2,-1,1,1,false,true,false)' > & ! Tbar_t-channel_lep.log 

echo 'DONE single top-bar t-channel'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2604,"2012_53x",2,-1,1,1,false,true,false)' > & ! T_tW-channel_lep.log 

echo 'DONE single top tW-channel'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2605,"2012_53x",2,-1,1,1,false,true,false)' > & ! Tbar_tW-channel_lep.log 

echo 'DONE single top-bar tW-channel'



echo 'Running W+jets samples....'

#root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2400,"2012_53x",2,-1,1,1,false,true,false)' > & ! WJetsToLNu_lep.log 

#echo 'DONE W+jets inclusive'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2401,"2012_53x",2,-1,1,1,false,true,false)' > & ! W1JetsToLNu_lep.log 

echo 'DONE W+1jet'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2402,"2012_53x",2,-1,1,1,false,true,false)' > & ! W2JetsToLNu_lep.log 

echo 'DONE W+2jets'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2403,"2012_53x",2,-1,1,1,false,true,false)' > & ! W3JetsToLNu_lep.log 

echo 'DONE W+3jets'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2404,"2012_53x",2,-1,1,1,false,true,false)' > & ! W4JetsToLNu_lep.log 

echo 'DONE W+4jets'



echo 'Running Z+jets samples....'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2850,"2012_53x",2,-1,1,1,false,true,false)' > & ! DYJetsToLL_M-10To50_lep.log 

echo 'DONE Z+jets M=10-50 inclusive'

#root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2800,"2012_53x",2,-1,1,1,false,true,false)' > & ! DYJetsToLL_M-50_lep.log 

#echo 'DONE Z+jets inclusive'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2801,"2012_53x",2,-1,1,1,false,true,false)' > & ! DY1JetsToLL_M-50_lep.log 

echo 'DONE Z+1jet'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2802,"2012_53x",2,-1,1,1,false,true,false)' > & ! DY2JetsToLL_M-50_lep.log 

echo 'DONE Z+2jets'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2803,"2012_53x",2,-1,1,1,false,true,false)' > & ! DY3JetsToLL_M-50_lep.log 

echo 'DONE Z+3jets'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2804,"2012_53x",2,-1,1,1,false,true,false)' > & ! DY4JetsToLL_M-50_lep.log 

echo 'DONE Z+4jets'



echo 'Running diboson samples....'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2700,"2012_53x",2,-1,1,1,false,true,false)' > & ! WW_lep.log 

echo 'DONE WW'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2701,"2012_53x",2,-1,1,1,false,true,false)' > & ! WZ_lep.log 

echo 'DONE WZ'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2702,"2012_53x",2,-1,1,1,false,true,false)' > & ! ZZ_lep.log 

echo 'DONE ZZ'



echo 'Running ttZ samples....'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2524,"2012_53x",2,-1,1,1,false,true,false)' > & ! TTWJets_lep.log 

echo 'DONE ttW+jets'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2523,"2012_53x",2,-1,1,1,false,true,false)' > & ! TTZJets_lep.log 

echo 'DONE ttZ+jets'



echo 'Running tt samples....'

#root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2500,"2012_53x",2,-1,1,1,false,true,false)' > & ! TTJets_MassiveBinDECAY_LF_lep.log 

#echo 'DONE tt+jets LF'

#root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2544,"2012_53x",2,-1,1,1,false,true,false)' > & ! TTJets_MassiveBinDECAY_CC_lep.log 

#echo 'DONE tt+jets CC'

#root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2555,"2012_53x",2,-1,1,1,false,true,false)' > & ! TTJets_MassiveBinDECAY_BB_lep.log 

#echo 'DONE tt+jets BB'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2566,"2012_53x",2,-1,1,1,false,true,false)' > & ! TTJets_HadronicMGDecays_LF_lep.log 

echo 'DONE hadronic tt+jets LF'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2576,"2012_53x",2,-1,1,1,false,true,false)' > & ! TTJets_HadronicMGDecays_CC_lep.log 

echo 'DONE hadronic tt+jets CC'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2586,"2012_53x",2,-1,1,1,false,true,false)' > & ! TTJets_HadronicMGDecays_BB_lep.log 

echo 'DONE hadronic tt+jets BB'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2563,"2012_53x",2,-1,1,1,false,true,false)' > & ! TTJets_SemiLeptMGDecays_LF_lep.log 

echo 'DONE semi-lep tt+jets LF'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2573,"2012_53x",2,-1,1,1,false,true,false)' > & ! TTJets_SemiLeptMGDecays_CC_lep.log 

echo 'DONE semi-lep tt+jets CC'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2583,"2012_53x",2,-1,1,1,false,true,false)' > & ! TTJets_SemiLeptMGDecays_BB_lep.log 

echo 'DONE semi-lep tt+jets BB'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2533,"2012_53x",2,-1,1,1,false,true,false)' > & ! TTJets_FullLeptMGDecays_LF_lep.log 

echo 'DONE full-lep tt+jets LF'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2543,"2012_53x",2,-1,1,1,false,true,false)' > & ! TTJets_FullLeptMGDecays_CC_lep.log 

echo 'DONE full-lep tt+jets CC'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(2553,"2012_53x",2,-1,1,1,false,true,false)' > & ! TTJets_FullLeptMGDecays_BB_lep.log 

echo 'DONE full-lep tt+jets BB'




echo 'Running data samples....'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(-1,"2012_53x",1,-1,1,1,false,true,false)' > & ! SingleMu_lep.log 

echo 'DONE single mu'

root -b -q head.C yggdrasil_treeReader_53x_v1.C+'(-2,"2012_53x",0,-1,1,1,false,true,false)' > & ! SingleEle_lep.log 

echo 'DONE single ele'

echo '............................DONE!'
