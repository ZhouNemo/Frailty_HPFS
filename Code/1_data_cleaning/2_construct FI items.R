# ==============================================================================
# Project: Frailty Trajectories Before and After Incident Cancer in the Health Professionals Follow-up Study
# Script: 2_construct FI items.R
# Author: Nemo Zhou
# Date started: Unknown (pre-existing script before documentation standard was applied)
# Date last updated: 2026-06-28
# Purpose: Constructs and scores individual frailty-index deficit items across questionnaire cycles before final renaming, imputation, and frailty-index calculation.
# ==============================================================================

library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(hpfs) # Internal package loader

# ==============================================================================
# 0. LOAD EXTERNAL DATASETS
# ==============================================================================
print("Loading Anthropometrics...")
anthropometrics_df <- load_hpfs_anthropometrics()

print("Loading Physical Activity...")
physact_df <- load_hpfs_physact()

# ==============================================================================
# 1. DEFINE PASSTHROUGH & VARIABLE MAPPINGS
# ==============================================================================

# --- A. Passthrough Variables (Smart NA) ---
# 1 = Participant skipped the entire section (Result should be NA)
disease_pt_map <- list(
  hp86 = NULL,       hp88 = "odxpt88",  hp90 = "q18pt90",  hp92 = "q23pt92", 
  hp94 = "q20pt94",  hp96 = "q25pt96",  hp98 = "q22pt98",  hp00 = "dxptb00", 
  hp02 = "dxptb02",  hp04 = "dxptb04",  hp06 = "dxptb06",  hp08 = "dxptb08", 
  hp10 = "dxptb10",  hp12 = "dxptb12",  hp14 = "dxptb14",  hp16 = "dxptb16", 
  hp18 = "dxptb18",  hp20 = "dxptb20"
)

med_pt_map <- list(
  hp86 = "q18pt86",  hp88 = "q26pt88",  hp90 = "q20pt90",  hp92 = NULL, 
  hp94 = "q22pt94",  hp96 = "q33pt96",  hp98 = "q23pt98",  hp00 = "medpt00", 
  hp02 = "medpt02",  hp04 = "medpt04",  hp06 = "medpt06",  hp08 = "medpt08", 
  hp10 = "medpt10",  hp12 = "medpt12",  hp14 = "medpt14",  hp16 = "medspt16", 
  hp18 = "medpt18",  hp20 = "medpt20"
)

# --- B. Disease Mappings ---
disease_maps <- list(
  cancer = list(
    hp86 = c("colc86", "lymp86", "mel86", "pros86", "ocan86"),
    hp88 = c("colc88", "lymp88", "mel88", "pros88", "ocan88"),
    hp90 = c("colc90", "lymp90", "mel90", "pros90", "ocan90"),
    hp92 = c("colc92", "lymp92", "mel92", "pros92", "ocan92"),
    hp94 = c("colc94", "lymp94", "mel94", "pros94", "ocan94"),
    hp96 = c("colc96", "lymp96", "mel96", "pros96", "ocan96"),
    hp98 = c("colc98", "lymp98", "mel98", "pros98", "ocan98"),
    hp00 = c("colc00", "lymp00", "mel00", "pros00", "ocan00"),
    hp02 = c("colc02", "lymp02", "mel02", "pros02", "ocan02"),
    hp04 = c("colc04", "lymp04", "mel04", "pros04", "ocan04"),
    hp06 = c("colc06", "lymp06", "mel06", "pros06", "ocan06"),
    hp08 = c("colc08", "lymp08", "mel08", "pros08", "ocan08", "blad08"),
    hp10 = c("colc10", "lymp10", "mel10", "pros10", "ocan10", "blad10"),
    hp12 = c("colc12", "lymp12", "mel12", "pros12", "ocan12", "blad12"),
    hp14 = c("colc14", "lymp14", "mel14", "pros14", "ocan14", "blad14"),
    hp16 = c("colc16", "lymp16", "mel16", "pros16", "ocan16", "blad16", "kidca16"),
    hp18 = c("colc18", "lymp18", "mel18", "pros18", "ocan18", "blad18", "kidca18"),
    hp20 = c("colc20", "lymp20", "mel20", "pros20", "ocan20", "blad20", "kidca20")
  ),
  chol = list(hp86="chol86", hp88="chol88", hp90="chol90", hp92="chol92", hp94="chol94", hp96="chol96", hp98="chol98", hp00="chol00", hp02="chol02", hp04="chol04", hp06="chol06", hp08="chol08", hp10="chol10", hp12="chol12", hp14="chol14", hp16="chol16", hp18="chol18", hp20="chol20"),
  hbp = list(hp86="hbp86", hp88="hbp88", hp90="hbp90", hp92="hbp92", hp94="hbp94", hp96="hbp96", hp98="hbp98", hp00="hbp00", hp02="hbp02", hp04="hbp04", hp06="hbp06", hp08="hbp08", hp10="hbp10", hp12="hbp12", hp14="hbp14", hp16="hbp16", hp18="hbp18", hp20="hbp20"),
  clau = list(hp86="clau86", hp88="clau88", hp90="clau90", hp92="clau92", hp94="clau94", hp96="clau96", hp98="clau98", hp00="clau00", hp02="clau02", hp04="clau04", hp06="clau06", hp08="clau08", hp10="clau10", hp12="clau12", hp14="clau14", hp16="peri16", hp18="peri18", hp20="peri20"),
  pe = list(hp86="pe86", hp88="pe88", hp90="pe90", hp92="pe92", hp94="pe94", hp96="pe96", hp98="pe98", hp00="pe00", hp02="pe02", hp04="pe04", hp06="pe06", hp08="pe08", hp10="pe10", hp12="pe12", hp14="pe14", hp16="pe16", hp18="pe18", hp20="pulm20"),
  visual = list(hp86="glau86", hp88=c("glau88", "macu88"), hp90=c("glau90", "macu90"), hp92=c("glau92", "macu92"), hp94=c("glau94", "macu94"), hp96=c("glau96", "macu96"), hp98=c("glau98", "macu98"), hp00=c("glau00", "macu00"), hp02=c("glau02", "macu02"), hp04=c("glau04", "macu04"), hp06=c("glau06", "macu06"), hp08=c("glau08", "macu08"), hp10=c("glau10", "macu10"), hp12=c("glau12", "macu12"), hp14=c("glau14", "macu14"), hp16=c("glau16", "macu16"), hp18=c("glau18", "macu18"), hp20=c("glau20", "macu20")),
  cad = list(hp86=c("ang86", "cabg86", "mi86"), hp88=c("ang88", "cabg88", "mi88"), hp90=c("ang90", "cabg90", "mi90"), hp92=c("ang92", "cabg92", "mi92"), hp94=c("ang94", "cabg94", "mi94"), hp96=c("ang96", "cabg96", "mi96"), hp98=c("ang98", "cabg98", "mi98"), hp00=c("ang00", "cabg00", "mi00"), hp02=c("ang02", "cabg02", "mi02"), hp04=c("ang04", "cabg04", "mi04"), hp06=c("ang06", "cabg06", "mi06"), hp08=c("ang08", "cabg08", "mi08"), hp10=c("ang10", "cabg10", "mi10"), hp12=c("ang12", "cabg12", "mi12"), hp14=c("ang14", "cabg14", "mi14"), hp16=c("ang16", "cabg16", "mi16"), hp18=c("ang18", "cabg18", "mi18"), hp20=c("ang20", "cabg20", "mi20")),
  cvd = list(hp86="str86", hp88="str88", hp90=c("str90", "cart90"), hp92=c("str92", "tia92", "cart92"), hp94=c("str94", "tia94", "cart94"), hp96=c("str96", "tia96", "cart96"), hp98=c("str98", "tia98", "cart98"), hp00=c("strk00", "tia00", "cart00"), hp02=c("strk02", "tia02", "cart02"), hp04=c("strk04", "tia04", "cart04"), hp06=c("strk06", "tia06", "cart06"), hp08=c("strk08", "tia08", "cart08"), hp10=c("strk10", "tia10", "cart10"), hp12=c("strk12", "tia12", "cart12"), hp14=c("strk14", "tia14", "cart14"), hp16=c("stk16", "tia16", "cart16"), hp18=c("strk18", "tia18", "cart18"), hp20=c("strk20", "tia20", "cart20")),
  gout = list(hp86="gout86", hp88="gout88", hp90="gout90", hp92="gout92", hp94="gout94", hp96="gout96", hp98="gout98", hp00="gout00", hp02="gout02", hp04="gout04", hp06="gout06", hp08="gout08", hp10="gout10", hp12="gout12", hp14="gout14", hp16="gout16", hp18="gout18", hp20="gout20"),
  hfra = list(hp86="frac86", hp88="hfra88", hp90="hfra90", hp92="hfra92", hp94="hfra94", hp96="hfra96", hp98="hfra98", hp00="hfra00", hp02="hfra02", hp04="hfra04", hp06="hfra06", hp08="hfra08", hp10="hfra10", hp12="hipf12", hp14="hipf14", hp16="hipf16", hp18="hipf18", hp20="hipf20"),
  ost = list(hp86="ost86", hp88="ost88", hp90="ost90", hp92="ost92", hp94="ost94", hp96="ost96", hp98="ost98", hp00="ost00", hp02="ost02", hp04="ost04", hp06="ost06", hp08="ost08", hp10="ost10", hp12="ost12", hp14="ost14", hp16="ost16", hp18="ost18", hp20="ost20"),
  arth = list(hp86="oa86", hp88="oa88", hp90="oa90", hp92="oa92", hp94="oa94", hp96="oa96", hp98="oa98", hp00="oa00", hp02="oa02", hp04="oa04", hp06="oa06", hp08="oa08", hp10="oa10", hp12="oa12", hp14="oa14", hp16="arth16", hp18="arth18", hp20="arth20"),
  dvrt = list(hp86=character(0), hp88=character(0), hp90="dvrt90", hp92="dvrt92", hp94="dvrt94", hp96="dvrt96", hp98="dvrt98", hp00="dvrt00", hp02="dvrt02", hp04="dvrt04", hp06="dvrt06", hp08="dvrt08", hp10="dvrt10", hp12="dvrt12", hp14="dvrt14", hp16="dvrt16", hp18="diver18", hp20="diver20"),
  ulcer = list(hp86="gast86", hp88="ulc88", hp90="gdul90", hp92="gdul92", hp94="gdul94", hp96="gdul96", hp98="gdul98", hp00="gdul00", hp02="gdul02", hp04="gdul04", hp06="gdul06", hp08="gdul08", hp10="gdul10", hp12="gdul12", hp14="gdul14", hp16="gdul16", hp18="gdul18", hp20="gdul20"),
  ucol = list(hp86="ucol86", hp88="ucol88", hp90="ucol90", hp92="ucol92", hp94="ucol94", hp96="ucol96", hp98="ucol98", hp00="ucol00", hp02="ucol02", hp04="ucol04", hp06="ucol06", hp08="ucol08", hp10="ucol10", hp12="ucol12", hp14="ucol14", hp16="ucol16", hp18="ucol18", hp20="ucol20"),
  park = list(hp86=character(0), hp88="park88", hp90="park90", hp92="park92", hp94="park94", hp96="park96", hp98="park98", hp00="park00", hp02="park02", hp04="park04", hp06="park06", hp08="park08", hp10="park10", hp12="park12", hp14="park14", hp16="park16", hp18="park18", hp20="park20"),
  db = list(hp86="db86", hp88="db88", hp90="db90", hp92="db92", hp94="db94", hp96="db96", hp98="db98", hp00="db00", hp02="db02", hp04="db04", hp06="db06", hp08="db08", hp10="db10", hp12="db12", hp14="db14", hp16="db16", hp18="db18", hp20="db20"),
  asthma = list(hp86="asth86", hp88=character(0), hp90=character(0), hp92="asth92", hp94=character(0), hp96="asth96", hp98="asth98", hp00="asth00", hp02="asth02", hp04="asth04", hp06="asth06", hp08="asth08", hp10="asth10", hp12="asth12", hp14="asth14", hp16=character(0), hp18="asth18", hp20=character(0))
)

# --- C. Medication Mappings (13 Classes) ---
poly_map <- list(
  Acetaminophen = list(hp86="tyl86", hp88="tyl88", hp90="tyl90", hp92="tyl92", hp94="tyl94", hp96="tyl96", hp98="tyl98", hp00="tyl00", hp02="tyl02", hp04="tyl04", hp06="tyl06", hp08="tyl08", hp10="tyl10", hp12="tyl12", hp14="tyl14", hp16="acetam16", hp18="tyl18", hp20="tyl20"),
  Aspirin = list(hp86="asp86", hp88="asp88", hp90="asp90", hp92="asp92", hp94="asp94", hp96="asp96", hp98="asp98", hp00="asp00", hp02="asp02", hp04="asp04", hp06="asp06", hp08="asp08", hp10="asp10", hp12="asp12", hp14=c("asp14","basp14"), hp16=c("asp16","baby16"), hp18=c("asp18","basp18"), hp20=c("asp20","basp20")),
  NSAIDs = list(hp86="motrn86", hp88="motrn88", hp90="motrn90", hp92="motrn92", hp94="motrn94", hp96=c("motrn96","onsai96"), hp98=c("motrn98","nsaid98"), hp00=c("motrn00","nsaid00"), hp02=c("mtrn02","cox2i02","analg02"), hp04=c("mtrn04","cox2i04","analg04"), hp06=c("mtrn06","cox2i06","analg06"), hp08=c("mtrn08","cox2i08","analg08"), hp10=c("mtrn10","cox2i10","analg10"), hp12=c("mtrn12","cox2i12","analg12"), hp14=c("mtrn14","cox2i14","analg14"), hp16=c("ibupro16","celebrex16","analg16"), hp18=c("mtrn18","cox2i18","analg18"), hp20=c("mtrn20","cox2i20","analg20")),
  BetaBlockers = list(hp86="betab86", hp88="betab88", hp90="betab90", hp92="betab92", hp94="betab94", hp96="betab96", hp98="betab98", hp00="betab00", hp02="betab02", hp04="betab04", hp06="betab06", hp08="betab08", hp10="betab10", hp12="betab12", hp14="betab14", hp16="bblock16", hp18="bblock18", hp20="bblock20"),
  CalciumBlockers = list(hp86="calcb86", hp88="ccblo88", hp90="ccblo90", hp92="ccblo92", hp94="ccblo94", hp96="ccblo96", hp98="ccblo98", hp00="ccblo00", hp02="calcb02", hp04="calcb04", hp06="calcb06", hp08="calcb08", hp10="calcb10", hp12="calcb12", hp14="calcb14", hp16="cblock16", hp18="cblock18", hp20="cblock20"),
  Diuretics = list(hp86=c("lasix86","diur86"), hp88=c("lasix88","diur88"), hp90=c("lasix90","diur90"), hp92=c("lasix92","thiaz92"), hp94=c("lasix94","thiaz94"), hp96=c("lasix96","thiaz96"), hp98=c("lasix98","thiaz98"), hp00=c("lasix00","thiaz00"), hp02=c("lasix02","thiaz02"), hp04=c("lasix04","thiaz04"), hp06=c("lasix06","thiaz06"), hp08=c("lasix08","thiaz08"), hp10=c("lasix10","thiaz10"), hp12=c("lasix12","thiaz12"), hp14=c("lasix14","thiaz14"), hp16=c("lasix16","thiaz16"), hp18=c("lasix18","thiaz18"), hp20="thiaz20"),
  OtherAntiHTN = list(hp86="ald86", hp88="ald88", hp90="ald90", hp92="ald92", hp94="ald94", hp96="ald96", hp98="ald98", hp00="oanth00", hp02="anthp02", hp04=c("ace04","anthp04"), hp06=c("ace06","anthp06"), hp08=c("ace08","arb08","anthp08"), hp10=c("ace10","arb10","anthp10"), hp12=c("ace12","arb12","anthp12"), hp14=c("ace14","arb14","anthp14"), hp16=c("aceinhb16","angio16","antihy16"), hp18=c("aceinhb18","angio18","antihy18"), hp20=c("aceinhb20","angio20","antihy20")),
  Cholesterol = list(hp86="antch86", hp88="chrx88", hp90="chrx90", hp92="chrx92", hp94="chrx94", hp96="chrx96", hp98="chrx98", hp00=c("chrx00","ochrx00"), hp02=c("stat02","ochrx02"), hp04=c("stat04","mev04","zoc04","crest04","prav04","lip04","ostat04","ochrx04"), hp06=c("stat06","mev06","prav06","zoc06","lip06","crest06","ostat06","ochrx06"), hp08=c("stat08","mev08","zoc08","crest08","prav08","lip08","ostat08","statpt08","ochrx08"), hp10=c("stat10","mev10","zoc10","crest10","prav10","lip10","ostat10","statpt10","ochrx10"), hp12=c("stat12","mev12","zoc12","crest12","prav12","lip12","ostat12","statpt12","ochrx12"), hp14=c("stat14","mev14","prav14","zoc14","lip14","crest14","ostat14","statpt14","ochrx14"), hp16=c("stat16","lipit16","prava16","crestor16","zocor16","othmed16","ochrx16"), hp18=c("stat18","lipit18","prava18","crestor18","zocor18","othmed18","ochrx18"), hp20=c("stat20","lipit20","prava20","crestor20","zocor20","othmed20","ochrx20")),
  AcidSupp = list(hp86="tag86", hp88="tag88", hp90="tag90", hp92="tag92", hp94="tag94", hp96="tag96", hp98="tag98", hp00="tag00", hp02="tag02", hp04=c("tag04","pril04"), hp06=c("tag06","pril06"), hp08=c("tag08","pril08"), hp10=c("tag10","pril10"), hp12=c("tag12","pril12"), hp14=c("tag14","pril14"), hp16=c("h2block16","prilo16"), hp18=c("h2block18","prilo18"), hp20=c("h2block20","prilo20")),
  Steroids = list(hp86=character(0), hp88="ster88", hp90="ster90", hp92="ster92", hp94="ster94", hp96="ster96", hp98="ster98", hp00="ster00", hp02="ster02", hp04="ster04", hp06="ster06", hp08="ster08", hp10="ster10", hp12="ster12", hp14="ster14", hp16="steroid16", hp18="steroid18", hp20="steroid20"),
  Antidepressants = list(hp86=character(0), hp88=character(0), hp90="antid90", hp92="antid92", hp94="antid94", hp96=c("przc96","tcyc96","antid96"), hp98=c("przc98","tcyc98","antid98"), hp00=c("przc00","tcyc00","antid00"), hp02=c("przc02","tcyc02","antid02"), hp04=c("przc04","tcyc04","antid04"), hp06=c("przc06","tcyc06","antid06"), hp08=c("ssri08","snri08","tcyc08","maoi08","antid08"), hp10="andep10", hp12=c("ssri12","andep12"), hp14=c("ssri14","andep14"), hp16=c("ssris16","tric16","antidep16"), hp18=c("ssris18","tric18","antidep18"), hp20=c("ssris20","tric20","antidep20")),
  Tranquilizers = list(hp86=character(0), hp88=character(0), hp90="tranq90", hp92=c("val92","thor92"), hp94=c("val94","thor94"), hp96=c("val96","thor96"), hp98="tranq98", hp00="val00", hp02="val02", hp04="val04", hp06="val06", hp08=c("benzo08","antpsy08"), hp10=character(0), hp12="tranq12", hp14="tranq14", hp16="minor_tranq16", hp18="minor_tranq18", hp20="minor_tranq20"),
  OtherRx = list(hp86=c("antar86","nitr86","digox86"), hp88=c("orx88","theo88","levo88","digox88","nitr88","antar88"), hp90=c("orx90","theo90","levo90","nitr90","digox90","antar90"), hp92=c("orx92","theo92","levo92","digox92","nitr92","antar92"), hp94="antar94", hp96=c("orx96","prosc96","alphb96","couma96","digox96"), hp98=c("orx98","prosc98","alphb98","couma98","digox98","omedp98"), hp00=c("orx00","prosc00","alphb00","couma00","digox00"), hp02=c("orx02","prosc02","alphb02","couma02","digox02"), hp04=c("orx04","prosc04","alphb04","couma04","digox04"), hp06=c("orx06","prosc06","alphb06","couma06","digox06"), hp08=c("orx08","antcon08","prosc08","alphb08","fosmx08","slpmed08","couma08","clopi08","insul08","ohypo08"), hp10=c("orx10","prosc10","alphb10","fosmx10","clopi10","couma10","slpmed10","insul10","ohypo10"), hp12=c("prosc12","prope12","avoda12","alphb12","fosmx12","couma12","clopi12","prada12","insul12","metfo12","avand12","hypog12"), hp14=c("antiarr14","prada14","insul14","metfo14","avand14","hypog14","prosc14","prope14","avoda14","alphb14","fosmx14","clopi14","couma14"), hp16=c("othermed16","aricept16","namenda16","proscar16","propecia16","avodart16","ablock16","ambien16","trazadone16","coum16","pradaxa16","clop16","digoxin16","antiarr16","insulin16","metformin16","actos16","hypogly16","opioid16"), hp18=c("othermed18","aricept18","namenda18","fosamax18","beta_agonists18","albuterol18","other_ba18","ambien18","trazadone18","coum18","pradaxa18","clop18","digoxin18","antiarr18","insulin18","metformin18","actos18","hypogly18","opioid18","proscar18","propecia18","avodart18","ablock18"), hp20=c("othermed20","aricept20","exelon20","razadyne20","namenda20","proscar20","propecia20","avodart20","ablock20","fosamax20","beta_agonists20","ambien20","otherslpmed20","coum20","pradaxa20","clop20","digoxin20","antiarr20","insulin20","ninsulinj20","metformin20","jardiance20","invokana20","farxiga20","januvia20","hypogly20","opioid20"))
)

# --- D. Functional & TV Mappings ---
func_maps <- list(
  stairs = list(
    hp86 = list(var=character(0), map=c()),
    hp88 = list(var="dact88", map=c("1"=0.0, "2"=1.0, "3"=NA)), 
    hp90 = list(var="dflt90", map=c("1"=0.0, "2"=1.0, "3"=NA)), 
    hp92 = list(var="dact92", map=c("1"=0.0, "2"=1.0, "3"=NA)), 
    hp94 = list(var="dact94", map=c("1"=0.0, "2"=1.0, "3"=NA)),
    hp96 = list(var="dflt96", map=c("1"=1.0, "2"=1.0, "3"=0.0, "4"=NA)), 
    hp98 = list(var="dact98", map=c("1"=0.0, "2"=1.0, "3"=NA)), 
    hp00 = list(var="dact00", map=c("1"=0.0, "2"=1.0, "3"=NA)), 
    hp02 = list(var="dact02", map=c("1"=0.0, "2"=1.0, "3"=NA)), 
    hp04 = list(var="dact04", map=c("1"=0.0, "2"=1.0, "3"=NA)),
    hp06 = list(var="dact06", map=c("1"=0.0, "2"=1.0, "3"=NA)), 
    hp08 = list(var="dact08", map=c("1"=0.0, "2"=1.0, "3"=NA)),
    hp10 = list(var=character(0), map=c()), 
    hp12 = list(var="diffst12", map=c("1"=0.0, "2"=1.0, "3"=NA)),
    hp14 = list(var="diffst14", map=c("1"=0.0, "2"=1.0, "3"=NA)), 
    hp16 = list(var="limonestair16", map=c("1"=1.0, "2"=1.0, "3"=0.0, "4"=NA)), 
    hp18 = list(var=character(0), map=c()), 
    hp20 = list(var="diffst20", map=c("1"=0.0, "2"=1.0, "3"=NA))
  ),
  balance = list(
    hp86 = list(var=character(0), map=c()), hp88 = list(var=character(0), map=c()),
    hp90 = list(var="dbal90", map=c("1"=0.0, "2"=1.0, "3"=NA)),
    hp92 = list(var="balnc92", map=c("1"=0.0, "2"=1.0, "3"=NA)),
    hp94 = list(var="balnc94", map=c("1"=0.0, "2"=1.0, "3"=NA)),
    hp96 = list(var="balnc96", map=c("1"=0.0, "2"=1.0, "3"=NA)),
    hp98 = list(var="balnc98", map=c("1"=0.0, "2"=1.0, "3"=NA)),
    hp00 = list(var="balnc00", map=c("1"=0.0, "2"=1.0, "3"=NA)),
    hp02 = list(var="balnc02", map=c("1"=0.0, "2"=1.0, "3"=NA)),
    hp04 = list(var="balnc04", map=c("1"=0.0, "2"=1.0, "3"=NA)),
    hp06 = list(var="balnc06", map=c("1"=0.0, "2"=1.0, "3"=NA)),
    hp08 = list(var="balnc08", map=c("1"=0.0, "2"=1.0, "3"=NA)),
    hp10 = list(var=character(0), map=c()),
    hp12 = list(var="balnc12", map=c("1"=0.0, "2"=1.0, "3"=NA)),
    hp14 = list(var=character(0), map=c()),
    hp16 = list(var="balnc16", map=c("1"=0.0, "2"=1.0, "3"=1.0, "4"=NA)),
    hp18 = list(var=character(0), map=c()),
    hp20 = list(var="balnc20", map=c("1"=0.0, "2"=1.0, "3"=1.0, "4"=NA))
  ),
  pace = list(
    hp86 = list(var="pace86", map=c("1"=1.0, "2"=0.0, "3"=0.0, "4"=0.0, "9"=NA)), 
    hp88 = list(var="pace88", map=c("1"=1.0, "2"=0.0, "3"=0.0, "4"=0.0, "9"=NA)), 
    hp90 = list(var=character(0), map=c()),
    hp92 = list(var="pace92", map=c("1"=1.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=NA)),
    hp94 = list(var="pace94", map=c("1"=1.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=NA)),
    hp96 = list(var="pace96", map=c("1"=1.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=NA)),
    hp98 = list(var="pace98", map=c("1"=1.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=NA)),
    hp00 = list(var="pace00", map=c("1"=1.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=NA)),
    hp02 = list(var=character(0), map=c()), hp04 = list(var=character(0), map=c()),
    hp06 = list(var=character(0), map=c()), hp08 = list(var=character(0), map=c()),
    hp10 = list(var="wpace10", map=c("1"=1.0, "2"=1.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=NA)), 
    hp12 = list(var="wpace12", map=c("1"=1.0, "2"=1.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=NA)), 
    hp14 = list(var="wpace14", map=c("1"=1.0, "2"=1.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=NA)), 
    hp16 = list(var="wpace16", map=c("1"=1.0, "2"=1.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=NA)),
    hp18 = list(var=character(0), map=c()), 
    hp20 = list(var="wpace20", map=c("1"=1.0, "2"=1.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=NA))
  )
)

tv_maps <- list(
  hp86 = c(),
  hp88 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.5, "6"=1.0, "9"=NA),
  hp90 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=0.0, "7"=0.0, "8"=0.0, "9"=0.0, "10"=0.0, "11"=0.5, "12"=0.5, "13"=1.0, "14"=NA, "99"=NA),
  hp92 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=0.0, "7"=0.0, "8"=0.0, "9"=0.0, "10"=0.0, "11"=0.5, "12"=0.5, "13"=1.0, "14"=NA, "99"=NA),
  hp94 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=0.0, "7"=0.0, "8"=0.0, "9"=0.0, "10"=0.0, "11"=0.5, "12"=0.5, "13"=1.0, "14"=NA, "99"=NA),
  hp96 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=0.0, "7"=0.0, "8"=0.0, "9"=0.0, "10"=0.0, "11"=0.5, "12"=0.5, "13"=1.0, "14"=NA, "99"=NA),
  hp98 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=0.0, "7"=0.0, "8"=0.0, "9"=0.0, "10"=0.0, "11"=0.5, "12"=0.5, "13"=1.0, "14"=NA, "99"=NA),
  hp00 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=0.0, "7"=0.0, "8"=0.0, "9"=0.0, "10"=0.0, "11"=0.5, "12"=0.5, "13"=1.0, "14"=NA, "99"=NA),
  hp02 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=0.0, "7"=0.0, "8"=0.0, "9"=0.0, "10"=0.0, "11"=0.5, "12"=0.5, "13"=1.0, "14"=NA, "99"=NA),
  hp04 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=0.0, "7"=0.0, "8"=0.0, "9"=0.0, "10"=0.0, "11"=0.5, "12"=0.5, "13"=1.0, "14"=NA, "99"=NA),
  hp06 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=0.0, "7"=0.0, "8"=0.0, "9"=0.0, "10"=0.0, "11"=0.5, "12"=0.5, "13"=1.0, "14"=NA, "99"=NA),
  hp08 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=0.0, "7"=0.0, "8"=0.0, "9"=0.0, "10"=0.0, "11"=0.5, "12"=0.5, "13"=1.0, "14"=NA, "99"=NA),
  hp10 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=0.0, "7"=0.0, "8"=0.0, "9"=0.0, "10"=0.0, "11"=0.5, "12"=0.5, "13"=1.0, "14"=NA, "99"=NA),
  hp12 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=0.0, "7"=0.0, "8"=0.0, "9"=0.0, "10"=0.0, "11"=0.5, "12"=0.5, "13"=1.0, "14"=NA, "99"=NA),
  hp14 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=0.0, "7"=0.0, "8"=0.0, "9"=0.0, "10"=0.0, "11"=0.5, "12"=0.5, "13"=1.0, "14"=NA, "99"=NA),
  hp16 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=0.5, "7"=1.0, "8"=1.0, "9"=1.0, "10"=NA),
  hp18 = c(),
  hp20 = c("1"=0.0, "2"=0.0, "3"=0.0, "4"=0.0, "5"=0.0, "6"=0.5, "7"=1.0, "8"=1.0, "9"=1.0, "10"=NA)
)

tv_vars <- list(hp86=character(0), hp88="sittv88", hp90="sittv90", hp92="sittv92", hp94="sittv94", hp96="sittv96", hp98="sittv98", hp00="sittv00", hp02="sittv02", hp04="sittv04", hp06="sittv06", hp08="sittv08", hp10="sittv10", hp12="sittv12", hp14="sittv14", hp16="sittv16", hp18=character(0), hp20="sittv20")


# ==============================================================================
# 2. MAIN EXTRACTION LOOP
# ==============================================================================
extracted_dfs <- list()
all_cycles <- paste0("hp", c("86", "88", "90", "92", "94", "96", "98", "00", "02", "04", "06", "08", "10", "12", "14", "16", "18", "20"))

for (cycle in all_cycles) {
  
  if(is.null(hp_data_list[[cycle]])) next
  df <- hp_data_list[[cycle]]
  yr <- substr(cycle, 3, 4)
  keep_cols <- c("id") 
  
  # ----------------------------------------------------------------------
  # PATCH 1: DERIVE MISSING ASPIRIN VAR (1996 and 1998)
  # ----------------------------------------------------------------------
  if (cycle == "hp96") {
    if ("aspd96" %in% names(df)) {
      df$asp96 <- case_when(
        df$aspd96 == 6 ~ NA_real_, df$aspd96 == 1 ~ 0.0, df$aspd96 %in% 2:5 ~ 1.0, TRUE ~ NA_real_
      )
    } else if ("nasp96" %in% names(df)) {
      df$asp96 <- case_when(
        df$nasp96 == 8 ~ NA_real_, df$nasp96 == 1 ~ 0.0, df$nasp96 %in% 2:7 ~ 1.0, TRUE ~ NA_real_
      )
    }
  }
  
  if (cycle == "hp98" && "aspfr98" %in% names(df)) {
    df$asp98 <- case_when(
      df$aspfr98 == 7 ~ NA_real_, 
      df$aspfr98 == 1 ~ 0.0, 
      df$aspfr98 %in% 2:6 ~ 1.0, 
      TRUE ~ NA_real_
    )
  }
  
  # ----------------------------------------------------------------------
  # PATCH 2: FIX DIVERTICULITIS CODING IN 2018 AND 2020 (1=No, 2=Yes)
  # ----------------------------------------------------------------------
  if (cycle %in% c("hp18", "hp20")) {
    diver_var <- paste0("diver", yr)
    if (diver_var %in% names(df)) {
      df[[diver_var]] <- case_when(
        df[[diver_var]] == 2 ~ 1.0, # Yes
        df[[diver_var]] == 1 ~ 0.0, # No
        TRUE ~ NA_real_
      )
    }
  }
  
  # --- 2A. PROCESS DISEASES WITH PASSTHROUGH LOGIC ---
  pt_var <- disease_pt_map[[cycle]]
  
  if (!is.null(pt_var) && pt_var %in% names(df)) {
    keep_cols <- c(keep_cols, pt_var)
    is_passthrough <- df[[pt_var]] == 1
  } else {
    is_passthrough <- rep(FALSE, nrow(df))
  }
  is_passthrough[is.na(is_passthrough)] <- FALSE
  
  for (c_type in names(disease_maps)) {
    vars_present <- intersect(disease_maps[[c_type]][[cycle]], names(df))
    keep_cols <- c(keep_cols, vars_present)
    
    new_var_name <- paste0(c_type, "_", yr)
    
    if (length(vars_present) > 0) {
      has_disease <- rowSums(df[vars_present] == 1, na.rm = TRUE) > 0
      
      df[[new_var_name]] <- case_when(
        is_passthrough ~ NA_real_,
        has_disease ~ 1.0,
        TRUE ~ 0.0    # <-- Missingness Logic: Blanks without passthrough evaluate to 0
      )
    } else {
      df[[new_var_name]] <- NA_real_
    }
  }
  
  # --- 2B. PROCESS MEDICATIONS WITH PASSTHROUGH LOGIC ---
  med_pt_var <- med_pt_map[[cycle]]
  
  if (!is.null(med_pt_var) && med_pt_var %in% names(df)) {
    keep_cols <- c(keep_cols, med_pt_var)
    is_med_passthrough <- df[[med_pt_var]] == 1
  } else {
    is_med_passthrough <- rep(FALSE, nrow(df))
  }
  is_med_passthrough[is.na(is_med_passthrough)] <- FALSE
  
  poly_flags <- list()
  
  for (cat_name in names(poly_map)) {
    cat_vars_present <- intersect(poly_map[[cat_name]][[cycle]], names(df))
    keep_cols <- c(keep_cols, cat_vars_present) 
    
    if (length(cat_vars_present) > 0) {
      has_med <- rowSums(df[cat_vars_present] == 1, na.rm = TRUE) > 0
      
      poly_flags[[cat_name]] <- case_when(
        is_med_passthrough ~ NA_integer_,
        has_med ~ 1L,
        TRUE ~ 0L      # <-- Missingness Logic: Blanks without passthrough evaluate to 0
      )
    } else {
      poly_flags[[cat_name]] <- NA_integer_
    }
  }
  
  df[[paste0("antidepressant_", yr)]] <- as.numeric(poly_flags[["Antidepressants"]])
  df[[paste0("tranq_", yr)]] <- as.numeric(poly_flags[["Tranquilizers"]])
  
  poly_matrix <- do.call(cbind, poly_flags)
  df[[paste0("polypharmacy_count_", yr)]] <- rowSums(poly_matrix, na.rm = TRUE)
  
  all_missing <- rowSums(is.na(poly_matrix)) == ncol(poly_matrix)
  df[[paste0("polypharmacy_count_", yr)]] <- if_else(all_missing, NA_real_, df[[paste0("polypharmacy_count_", yr)]])
  
  df[[paste0("polypharmacy_", yr)]] <- case_when(
    is.na(df[[paste0("polypharmacy_count_", yr)]]) ~ NA_real_,
    df[[paste0("polypharmacy_count_", yr)]] >= 5 ~ 1.0,
    TRUE ~ 0.0
  )
  
  # --- 2C. FUNCTIONAL CAPACITIES ---
  for (f_type in names(func_maps)) {
    var_info <- func_maps[[f_type]][[cycle]]
    if (length(var_info$var) > 0 && var_info$var %in% names(df)) {
      raw_col <- var_info$var
      map_vec <- var_info$map
      keep_cols <- c(keep_cols, raw_col)
      df[[paste0(f_type, "_", yr)]] <- map_vec[as.character(df[[raw_col]])]
    } else {
      df[[paste0(f_type, "_", yr)]] <- NA_real_
    }
  }
  
  # --- 2D. SITTING AT TV ---
  tv_raw_var <- tv_vars[[cycle]]
  if (length(tv_raw_var) > 0 && tv_raw_var %in% names(df)) {
    keep_cols <- c(keep_cols, tv_raw_var)
    current_tv_map <- tv_maps[[cycle]]
    df[[paste0("tv_score_", yr)]] <- current_tv_map[as.character(df[[tv_raw_var]])]
  } else {
    df[[paste0("tv_score_", yr)]] <- NA_real_
  }
  
  # --- Final Selection ---
  new_vars <- c(
    paste0(names(disease_maps), "_", yr),
    paste0(names(func_maps), "_", yr),
    paste0("antidepressant_", yr), 
    paste0("tranq_", yr), 
    paste0("polypharmacy_count_", yr), 
    paste0("polypharmacy_", yr),
    paste0("tv_score_", yr)
  )
  
  all_target_cols <- unique(c(keep_cols, new_vars))
  final_cols <- intersect(all_target_cols, names(df))
  extracted_dfs[[cycle]] <- df %>% select(all_of(final_cols))
}

# ==============================================================================
# 3. MERGE DATAFRAMES AND POST-PROCESS
# ==============================================================================
FI <- reduce(extracted_dfs, full_join, by = "id")

# Force FI id to character to guarantee a perfect match for future merges
FI <- FI %>% 
  select(-matches("^(act|act_score|bmi|abnormal_bmi_score|weight|wtloss_score)_")) %>%
  mutate(id = as.character(id))


# --- 3A. Process Physical Activity ---
act_wide <- physact_df %>%
  mutate(
    cycle_str = str_extract(as.character(cycle), "\\d{2}$"),
    id = as.character(id) 
  ) %>%
  filter(cycle_str %in% c("86", "88", "90", "92", "94", "96", "98", "00", "02", "04", "06", "08", "10", "12", "14", "16", "18", "20")) %>%
  mutate(
    act_score = case_when(
      is.na(act) ~ NA_real_,
      act < 3 ~ 1.0,
      act >= 3 & act <= 7.5 ~ 0.5,
      act > 7.5 ~ 0.0
    )
  ) %>%
  select(id, cycle_str, act, act_score) %>%
  pivot_wider(
    names_from = cycle_str, 
    values_from = c(act, act_score), 
    names_sep = "_" 
  )

# --- 3B. Process BMI & Weight ---
anthro_clean <- anthropometrics_df %>%
  mutate(
    cycle_str = str_extract(as.character(cycle), "\\d{2}$"),
    id = as.character(id)
  )

bmi_wide <- anthro_clean %>%
  filter(cycle_str %in% c("86", "88", "90", "92", "94", "96", "98", "00", "02", "04", "06", "08", "10", "12", "14", "16", "18", "20")) %>%
  mutate(
    abnormal_bmi_score = case_when(
      is.na(bmi) ~ NA_real_,
      bmi > 30 | bmi < 18.5 ~ 1.0,
      TRUE ~ 0.0
    )
  ) %>%
  select(id, cycle_str, bmi, abnormal_bmi_score) %>%
  pivot_wider(
    names_from = cycle_str, 
    values_from = c(bmi, abnormal_bmi_score), 
    names_sep = "_" 
  )

# Ensure baseline weight (84) is pulled if it exists for the 86 calculation
weight_wide <- anthro_clean %>%
  filter(cycle_str %in% c("84", "86", "88", "90", "92", "94", "96", "98", "00", "02", "04", "06", "08", "10", "12", "14", "16", "18", "20")) %>%
  select(id, cycle_str, wt) %>%
  pivot_wider(
    names_from = cycle_str, 
    values_from = wt, 
    names_prefix = "weight_" 
  )

# --- 3C. Final Merge ---
FI <- FI %>%
  left_join(act_wide, by = "id") %>%
  left_join(bmi_wide, by = "id") %>%
  left_join(weight_wide, by = "id")

# Ensure all required weight columns exist to prevent dplyr evaluation errors
expected_weights <- paste0("weight_", c("84", "86", "88", "90", "92", "94", "96", "98", "00", "02", "04", "06", "08", "10", "12", "14", "16", "18", "20"))

for (w in expected_weights) {
  if (!w %in% names(FI)) {
    FI[[w]] <- NA_real_ 
  }
}

# --- 3D. Calculate Weight Loss (>= 5%) Dynamically over 2-year periods ---
target_years <- c("86", "88", "90", "92", "94", "96", "98", "00", "02", "04", "06", "08", "10", "12", "14", "16", "18", "20")
prev_years   <- c("84", "86", "88", "90", "92", "94", "96", "98", "00", "02", "04", "06", "08", "10", "12", "14", "16", "18")

for (i in seq_along(target_years)) {
  curr <- target_years[i]
  prev <- prev_years[i]
  
  curr_wt <- paste0("weight_", curr)
  prev_wt <- paste0("weight_", prev)
  score_col <- paste0("wtloss_score_", curr)
  
  FI[[score_col]] <- case_when(
    is.na(FI[[prev_wt]]) | is.na(FI[[curr_wt]]) ~ NA_real_,
    FI[[prev_wt]] > 0 & ((FI[[prev_wt]] - FI[[curr_wt]]) / FI[[prev_wt]] >= 0.05) ~ 1.0,
    TRUE ~ 0.0
  )
}

print("Merge complete! Physical Activity, BMI, and Weight Loss successfully integrated.")

# ==============================================================================
# 4. SAVE OUTPUT
# ==============================================================================
target_dir <- "/n/home06/xyzhou/Frailty"

if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE)
}

write.csv(FI, file = file.path(target_dir, "FI_longitudinal_1986_2020.csv"), row.names = FALSE)
saveRDS(FI, file = file.path(target_dir, "FI_longitudinal_1986_2020.rds"))

cat("Files successfully saved to:", target_dir)