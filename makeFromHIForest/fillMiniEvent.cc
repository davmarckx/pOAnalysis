#include "interface/MiniEvent.h"
#include "interface/HIFEvent.h"
#include "interface/HIFEvent_gen.h"

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TDirectory.h"

#include "TH1.h"
#include "TH2D.h"
#include <cmath>
#include <iostream>

double deltaR(double eta1, double eta2, double phi1, double phi2){
  double deltaEta = eta1-eta2;
  double deltaPhi = phi1-phi2;
  if (deltaPhi > 3.1415){deltaPhi -= 6.2830;}
  if (deltaPhi < -3.1415){deltaPhi += 6.2830;}


  return std::sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
}

double fit_extended(double pval,double a,double b,double c,double d,double e,double f,double g,double h,double i){
    return  a/(pow(pval+g,4)) + e/(pow(pval+h,3))+ b/(pow(pval+i,2)) + c/(pval+f) + d;
}



// hardcoded fit results (dumped in dEdx_fits.txt, no longer used but kept in case)

std::pair<std::vector<double>,std::vector<double>> get_vminmax(std::string pID){
  std::vector<double> vmin;
  std::vector<double> vmax;

  if(pID == "pion"){
    vmin = {2.5522571182250977, 2.5522571182250977, 2.990761637687683, 2.8064851188659667, 4.413767585754394, 4.483224391937256, 3.680380392074585, 3.325924003124237, 3.081802132129669, 2.8715845346450806, 2.7359977436065672, 2.592973463535309, 2.5134663343429566, 2.442749423980713, 2.4063335990905763, 2.3530305099487303, 2.332056164741516, 2.2996073532104493, 2.277120132446289, 2.2657714819908144, 2.2539173793792724, 2.2489497423172, 2.2405522871017456, 2.2423375177383424, 2.234165072441101, 2.2400204658508303, 2.2295066499710083, 2.227503623962402, 2.2274788284301756, 2.230953822135925, 2.2353474521636962, 2.239715185165405, 2.238243145942688, 2.241450157165527, 2.243664562702179, 2.2414504861831666, 2.2506228399276735, 2.250705614089966, 2.2595553541183473, 2.258083658218384, 2.2589592218399046, 2.258362123966217, 2.265838556289673, 2.253099312782288, 2.271215534210205, 2.264350595474243, 2.2671510219573974, 2.270364480018616, 2.2672643280029297, 2.278645372390747, 2.2809073877334596, 2.278449251651764, 2.2823112487792967, 2.286814079284668, 2.2871182560920715, 2.2894587302207947, 2.2886719703674316, 2.291870903968811, 2.2846816253662108, 2.2927920770645143, 2.2949837493896483, 2.297338273525238, 2.3009276247024535, 2.297920460700989, 2.300873041152954, 2.3085290575027466, 2.2951488208770754, 2.30506112575531, 2.317436099052429, 2.316680612564087, 2.311418237686157, 2.3258571887016295, 2.3115813302993775, 2.3123785996437074, 2.319012999534607, 2.3247218990325926, 2.326101541519165, 2.3299714016914366, 2.321194922924042, 2.3320268392562866, 2.3237445068359377, 2.3294408917427063, 2.3340168952941895, 2.3416407394409178, 2.3355664730072023, 2.3297407269477843, 2.3442344665527344, 2.3457939434051513, 2.350414681434631, 2.338681221008301, 2.3370994424819944, 2.3462705206871033, 2.343661251068115, 2.336158494949341, 2.3403678560256957, 2.3630341696739197, 2.3507246589660644, 2.352370254993439, 2.3840819072723387};
    vmax = {27.707912139892578, 27.707912139892578, 27.707912139892578, 27.707912139892578, 19.709307479858403, 14.90684654235839, 11.236470928192134, 8.75263124465944, 7.139702453613275, 6.41166206359863, 5.735093555450439, 5.268970584869376, 4.94947519302368, 4.728464870452879, 4.538455905914307, 4.374221172332761, 4.174395799636841, 4.057551021575928, 3.9169489669799784, 3.85960078716278, 3.732796058654785, 3.67362494468689, 3.6255632305145262, 3.5968482875823984, 3.5351059436798096, 3.4982347869873056, 3.4578948116302493, 3.4220904874801645, 3.4306385231018024, 3.45078759670258, 3.396824111938477, 3.380784158706665, 3.404141216278075, 3.321674757003784, 3.3515714406967145, 3.3570011615753175, 3.345621814727783, 3.34949164390564, 3.3490824317932093, 3.3376774883270266, 3.3327192783355697, 3.2999357795715327, 3.350254201889038, 3.3270253753662105, 3.323412795066836, 3.310732202529905, 3.346562976837156, 3.393259501457213, 3.4144734573364257, 3.3638251781463615, 3.2909772491455076, 3.3653292274475075, 3.3017428112030025, 3.376480541229248, 3.399247646331787, 3.3359765768051113, 3.346496105194092, 3.3383659362792937, 3.33710391998291, 3.3456777763366685, 3.349086189270019, 3.3423323106765745, 3.2968100357055667, 3.320580692291259, 3.3083676052093494, 3.2804401826858545, 3.3177005386352545, 3.2758921861648553, 3.275201916694641, 3.2749911594390855, 3.245606679916382, 3.235558400154114, 3.217060966491698, 3.2207544040679927, 3.2358856201171875, 3.235369148254394, 3.2438872814178468, 3.2274667549133302, 3.1949121236801137, 3.197302598953247, 3.194299697875976, 3.1995123624801636, 3.188454771041869, 3.1859739303588857, 3.1817539691925036, 3.173678541183472, 3.169371280670166, 3.184587888717651, 3.1977388620376583, 3.160354971885681, 3.1720952987670907, 3.1357075452804564, 3.1502152252197266, 3.146205253601074, 3.1726896667480466, 3.1676759099960328, 3.1746576499938963, 3.187950377464294, 3.1848096752166746};
  }
  else if (pID == "kaon"){
    vmin = {49.0, 44.0, 39.0, 34.0, 29.0, 24.0, 20.0, 16.0, 14.0, 10.797355918884277, 10.325242462158203, 9.167609205245972, 8.403738632202149, 7.457221384048462, 6.758189878463745, 6.457902822494507, 5.693963956832886, 5.126575875282287, 4.88536235332489, 4.478871631622314, 4.2123742771148684, 4.088614706993103, 3.886824607849121, 3.798894205093384, 3.677508707046509, 3.570455479621887, 3.449812984466553, 3.299306819438934, 3.2692803859710695, 3.1726276206970216, 3.063761563301086, 3.052488408088684, 2.943612480163574, 2.9417372035980223, 2.8446310997009276, 2.8070102024078367, 2.7869019889831543, 2.7263815641403197, 2.713680930137634, 2.652940773963928, 2.6418435049057005, 2.641500473022461, 2.6103971505165102, 2.5938216257095337, 2.5849129915237428, 2.5590945053100587, 2.554396755695343, 2.5094133973121644, 2.5112310528755186, 2.516063942909241, 2.475096197128296, 2.4754973888397216, 2.482828698158264, 2.437875053882599, 2.448650588989258, 2.4539449620246887, 2.4470333528518675, 2.413846216201782, 2.413013582229614, 2.408626832962036, 2.4052226543426514, 2.4248095035552977, 2.398063292503357, 2.381012969017029, 2.3921443033218384, 2.372892417907715, 2.3618753242492674, 2.3475664162635805, 2.379637403488159, 2.3584392929077147, 2.3568777322769163, 2.33005713224411, 2.3453455567359924, 2.345520100593567, 2.3328378796577454, 2.3191423797607422, 2.332850389480591, 2.348707365989685, 2.3386594796180726, 2.3264032340049745, 2.335143766403198, 2.329082088470459, 2.313366425037384, 2.325324201583862, 2.3073085260391237, 2.3145986199378967, 2.2938318967819216, 2.324183406829834, 2.3050995445251465, 2.3137952709197998, 2.3306115746498106, 2.3323495388031006, 2.3259159421920774, 2.335581305027008, 2.3315239000320434, 2.2958042860031127, 2.304254298210144, 2.3310670232772828, 2.3064167141914367};
    vmax = {65.0, 60.0, 55.0, 50.0, 45.0, 40.0, 36.0, 32.0, 30.0, 26.42053665161133, 24.711823959350582, 21.44807540893553, 20.579441394805905, 17.878064689636233, 16.792983798980714, 13.767140827178956, 13.469088268280029, 11.86005916595459, 9.793175344467164, 8.573926067352298, 7.489092607498167, 7.153995456695557, 6.484981346130371, 6.025316867828369, 5.922265796661377, 5.521693921089171, 5.279190168380738, 5.067882108688354, 4.929754676818848, 4.770280513763429, 4.5744168663024904, 4.460778403282165, 4.250589141845703, 4.276648759841916, 4.141985416412352, 4.107042398452759, 3.999547433853149, 3.953354644775389, 3.9168520379066467, 3.8503212451934807, 3.818150327205658, 3.7534563541412354, 3.638160562515259, 3.6561514043807968, 3.60795158624649, 3.611114068031311, 3.596044912338257, 3.6050407886505127, 3.574023909568787, 3.4636739587783807, 3.3807236099243165, 3.506568193435668, 3.508215715885162, 3.388421697616577, 3.393072102069852, 3.3743178081512437, 3.2841577720642086, 3.3898258090019224, 3.339597969055175, 3.284492387771606, 3.229099988937378, 3.3036896324157716, 3.2539861345291134, 3.2152321743965144, 3.2430321645736693, 3.289692673683166, 3.252127389907837, 3.2111606645584105, 3.2255835056304925, 3.136999254226684, 3.2216280698776245, 3.17496684551239, 3.1518825769424432, 3.157311165332794, 3.167756676673889, 3.098926753997803, 3.1280198574066156, 3.10588104724884, 3.109976375102997, 3.1047714376449593, 3.105818061828613, 3.1196437454223624, 3.137750494480132, 3.1149710655212406, 3.1148919773101804, 3.069898424148559, 3.031081438064575, 3.067927207946777, 3.062179079055786, 3.0243753767013546, 3.0778401255607606, 3.0117238998413085, 3.068771717548371, 3.0986820650100704, 3.013539161682129, 3.075511372089385, 3.0576591682434078, 3.0338739871978766, 3.1222723484039308};
  }
  else if (pID == "proton"){
    vmin = {65.5, 60.5, 55.5, 53.5, 53.5, 53.0, 48.0, 43.0, 38.0, 33.0, 28.0, 24.0, 20.0, 18.0, 17.256237428477316, 13.0, 11.578316097259522, 10.696618347167968, 10.276598663330079, 9.9087912940979, 9.290714511871338, 9.0466255569458, 8.57047477722168, 8.081250953674317, 7.639736137390137, 7.36602599143982, 6.980116271972657, 6.710602450370788, 6.353794994354248, 6.1682517528533936, 5.843960485458374, 5.676335220336914, 5.534111738204956, 5.319490933418274, 5.167079148292541, 5.002064228057861, 4.878883390426636, 4.710498485565186, 4.575484037399292, 4.501122832298279, 4.358710737228393, 4.273690986633301, 4.135044975280762, 4.050842227935791, 3.957529091835022, 3.912015104293823, 3.8294858169555663, 3.78310170173645, 3.615698938369751, 3.5893679308891295, 3.557840690612793, 3.302420938014984, 3.326082239151001, 3.1420379281044006, 2.925914692878723, 2.78452299118042, 2.747127962112427, 2.693664507865906, 2.6762944078445434, 2.7005292534828187, 2.705940160751343, 2.699184477329254, 2.6390406489372253, 2.740687565803528, 2.749296314716339, 2.698706502914429, 2.6908410000801086, 2.6320170092582704, 2.6781671953201296, 2.6325659680366518, 2.635861396789551, 2.6629924464225767, 2.6174585747718813, 2.669738087654114, 2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727,2.686386227607727};
    vmax = {90.0,90.0,90.0,90.0,90.0,90.0,90.0,90.0,90.0,90.0,90.0,90.0,90.0,80.0,70.0, 60.0, 32.0, 30.1, 29.071761734008692, 25.67537796020506, 26.126737060546873, 24.528037109375017, 20.346788101196278, 19.78708869934082, 17.967738342285152, 17.33329048156739, 16.45245220184327, 15.577755470275878, 13.693615837097168, 13.060929775238037, 12.624228782653804, 11.116469459533665, 9.901309490203857, 9.222109413146985, 9.258628196716307, 7.782722043991086, 8.251549339294431, 7.5790978622436525, 7.035638093948365, 7.280358648300174, 6.795026702880858, 6.458363580703735, 6.0845237064361575, 6.074124050140381, 5.909504861831662, 5.727134132385251, 5.6849892902374295, 5.340808343887328, 5.431743965148925, 5.417793769836424, 5.028053817749024, 5.043228044509887, 4.969036979675293, 4.816554069519044, 4.828817377090454, 4.755254526138307, 4.525655822753905, 4.636913070678706, 4.52254383087158, 4.341978054046631, 4.3118659782409665, 4.180289030075073, 4.232295989990234, 4.259139099121095, 4.235232610702511, 4.195414686203, 4.074352703094482, 4.0283683395385745, 4.002374544143676, 3.9587430906295777, 4.132192325592041, 3.85005774974823, 3.928359794616698, 3.93667248249054, 3.7298201274871827, 3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827,3.7298201274871827};
  }
  else{
    std::cerr<<"pID not recognized in binned dE/dx fit.";
  }
  return std::make_pair(vmin,vmax);
}

// ------------ method called to find the max and min dE/dx values to be in the requested band  ------------

std::pair<double,double> get_minmax_binned(double pval,std::pair<std::vector<double>,std::vector<double>> pair){
 std::vector<double> vmin = pair.first;
 std::vector<double> vmax = pair.second;

 int nr = int(pval/0.0202);
 if (nr>98){
  nr = -1;
 }
 return std::make_pair(vmin[nr],vmax[nr]);
}

// ------------ method called to find the max and min dE/dx values to be in the requested band  ------------

std::pair<std::pair<double,double>,std::pair<double,double>> get_minmax(double pval,std::string flavor){
 double min;
 double max;
 double min_l;
 double max_l;
 if (flavor=="pion"){
    min = fit_extended(pval,0.40989992,  2.06147096, -1.54500484,  2.73634483, -1.41475266,  0.35678269,0.31638913,  0.31901972,  0.33499779);
    max = fit_extended(pval,4.48163187e+09, -1.99092710e+03,  8.28637374e+02, -5.04765327e+02, -1.75070459e+02,  2.03157757e-01, -6.17986544e+01,  4.33938800e-01,  1.10524518e+00);

    min = std::min(5.,min);
    max = std::min(27.5,max);

    min_l = fit_extended(pval+0.1,0.40989992,  2.06147096, -1.54500484,  2.73634483, -1.41475266,  0.35678269,0.31638913,  0.31901972,  0.33499779);
    max_l = fit_extended(pval-0.1,4.48163187e+09, -1.99092710e+03,  8.28637374e+02, -5.04765327e+02, -1.75070459e+02,  2.03157757e-01, -6.17986544e+01,  4.33938800e-01,  1.10524518e+00);

    min_l = std::min(4.,min_l);
    max_l = std::min(28.5,max_l);
  }
 else if (flavor=="kaon"){
    min = fit_extended(pval,24.23585851,  53.72213504, -27.85200971,   7.70967802, -48.64089385, 0.25341064,   0.32862305,   0.26308091,   0.26864415);
    max = fit_extended(pval,1.85534645e+03, -1.79032095e+06,  4.29904348e+04,  3.61034274e+03,  2.49962977e+02,  3.87791962e+00, -7.30777857e+01, -2.20433295e+00,  1.10700053e+01);

    min = std::min(40.,min);
    max = std::min(40.,max);

    min_l = fit_extended(pval+0.1,24.23585851,  53.72213504, -27.85200971,   7.70967802, -48.64089385, 0.25341064,   0.32862305,   0.26308091,   0.26864415);
    max_l = fit_extended(pval-0.1,1.85534645e+03, -1.79032095e+06,  4.29904348e+04,  3.61034274e+03,  2.49962977e+02,  3.87791962e+00, -7.30777857e+01, -2.20433295e+00,  1.10700053e+01);

    min_l = std::min(39.,min_l);
    max_l = std::min(41.,max_l);
  }
 else if (flavor=="proton"){
    double maxk = fit_extended(pval,1.85534645e+03, -1.79032095e+06,  4.29904348e+04,  3.61034274e+03,  2.49962977e+02,  3.87791962e+00, -7.30777857e+01, -2.20433295e+00,  1.10700053e+01);
    min = fit_extended(pval,9.19955745e+00,  8.92033952e+01, -5.63296835e+01,  1.72842046e+01,-3.46057114e+01, -7.06307739e-02,  1.56179518e-01,  1.61667844e-01,1.82998287e-01);
    max = fit_extended(pval,-2.40236732e+03, -1.08477323e+05,  2.49771514e+04, -2.63751597e+03,-4.83004264e+04,  1.76347258e+00 , 2.46524720e+03, -8.06810115e+00,3.07529625e+00);

    min = std::min(82.,min);
    min = std::max(maxk+0.1,min);
    max = std::min(150.,max);

    min_l = fit_extended(pval+0.1,9.19955745e+00,  8.92033952e+01, -5.63296835e+01,  1.72842046e+01,-3.46057114e+01, -7.06307739e-02,  1.56179518e-01,  1.61667844e-01,1.82998287e-01);
    max_l = fit_extended(pval-0.1,-2.40236732e+03, -1.08477323e+05,  2.49771514e+04, -2.63751597e+03,-4.83004264e+04,  1.76347258e+00 , 2.46524720e+03, -8.06810115e+00,3.07529625e+00);

    min_l = std::min(81.,min_l);
    max_l = std::min(151.,max_l);
  }
 else if (flavor=="deuterion"){
    double maxpr = fit_extended(pval,-2.40236732e+03, -1.08477323e+05,  2.49771514e+04, -2.63751597e+03,-4.83004264e+04,  1.76347258e+00 , 2.46524720e+03, -8.06810115e+00,3.07529625e+00);
    min = -0.1+fit_extended(pval,-5.18251880e+04,  6.00143418e+04,  1.89169865e+03, -6.36742783e+02,1.53154564e+03,  3.36913031e+00,  4.23613847e+00,  1.92057065e+02,-1.56391160e+01);
    max = 0.1+fit_extended(pval,1.80313156e+03,  6.67481452e+03, -2.94162671e+03,  2.61371955e+02, 6.01056572e+03,  4.28953682e+00,  8.35549262e+00,  6.95408925e+00, 3.73133240e+00);

    min = std::min(82.,min);
    min = std::max(maxpr+0.1,min);
    max = std::min(150.,max);

    min_l = -0.1+fit_extended(pval+0.1, -5.18251880e+04,  5.98043418e+04,  1.89169865e+03, -6.36742783e+02-0.1,1.53154564e+03,  3.36913031e+00,  4.23613847e+00,  1.92057065e+02,-1.56391160e+01);
    max_l = 0.1+fit_extended(pval-0.1, 1.80313156e+03,  6.67481452e+03, -2.94162671e+03,  2.61371955e+02+1.5, 6.01056572e+03,  4.28953682e+00,  8.35549262e+00,  6.95408925e+00, 3.73133240e+00);

    min_l = std::min(81.,min_l);
    min_l = std::max(maxpr,min_l);
    max_l = std::min(151.,max_l);
  }
  return  std::make_pair(std::make_pair(min,max),std::make_pair(min_l,max_l));
}



int main( int argc, char* argv[] ){

    std::cerr << "###starting###" << std::endl;

    int nargs = 2;
    if( argc != nargs+1 ){
        std::cerr << "ERROR: fillMiniEvent.cc requires " << std::to_string(nargs);
        std::cerr << " arguments to run...: " << std::endl;
        std::cerr << "  - inputfile" << std::endl;
        std::cerr << "  - outputname" << std::endl;
        return -1;
    }

    // parse arguments
    std::vector< std::string > argvStr( &argv[0], &argv[0] + argc );
    std::string& inputfile = argvStr[1];
    std::string& outputname = argvStr[2];

    // print arguments
    std::cout << "Found following arguments:" << std::endl;
    std::cout << "  - input file: " << inputfile << std::endl;
    std::cout << "  - output name (e.g. eposlhcr8): " << outputname << std::endl;


    // make output tree
    //std::string outputfile = "/eos/user/d/dmarckx/PO/CMSSW_15_0_0_pre2/src/pOAnalysis/OUTPUT/output_" + outputname +".root";
    std::string outputfile = "/pnfs/iihe/cms/store/user/dmarckx/pO_miniAOD/MiniEvent/output_" + outputname +"_withD.root";

    TFile of(outputfile.c_str(), "recreate");
    TDirectory *dir = of.mkdir("analysis");
    TTree *tree_ = new TTree("tree","tree");
    //tree_->SetDirectory(0);
    //tree_->SetName("tree");
    MiniEvent_t ev_;
    createMiniEventTree(tree_,ev_);
    

    std::cout<<"made input tree"<<std::endl;
    
    
    // make HIForest tree
    //TFile f(inputfile.c_str(), "read");

    std::cout<<"opening file "<< inputfile <<std::endl;


    HIFEvent_t hifev_;
    HIFEventgen_t hifgenev_;

    //TTree *HIFtree_ = f.Get<TTree>("PbPbTracks/trackTree");
    //TTree *HIFtreegen_ = f.Get<TTree>("HiGenParticleAna/hi");
    TChain *HIFtree_ = new TChain("PbPbTracks/trackTree");
    TChain *HIFtreegen_ = new TChain("HiGenParticleAna/hi");

    const std::string& recotree = "*";
    const std::string& gentree = "*";

    HIFtree_->Add((inputfile).c_str());
    HIFtreegen_->Add((inputfile).c_str());

    attachToHIFEventTree(HIFtree_,hifev_);
    attachToHIFEventTree_gen(HIFtreegen_,hifgenev_); 

    // init sample wide stuff
    ev_.isData  = false;
    ev_.run = hifev_.run;
    ev_.lumi = hifev_.lumi;

    // go over tracks
    ev_.ntrk = 0;
    Float_t dRmatch = 99.; // added
    Float_t pmatch = 99.; // added
    Int_t matchedPdgId = 0; // added
    bool hasReco = false; // added
    std::vector<int> v = {22,12,14,16,130,-12,-14,-16,-130}; // skip PiDs that do not leave tracks (photons, neutrinos, K0's)

    // init dEdx fit info
    std::pair<std::vector<double>,std::vector<double>> pionmap = get_vminmax("pion");
    std::pair<std::vector<double>,std::vector<double>> kaonmap = get_vminmax("kaon");
    std::pair<std::vector<double>,std::vector<double>> protonmap = get_vminmax("proton");

    // fit method
    std::pair<std::pair<double,double>,std::pair<double,double>> pionpair;
    std::pair<std::pair<double,double>,std::pair<double,double>> kaonpair;
    std::pair<std::pair<double,double>,std::pair<double,double>> protonpair;
    std::pair<std::pair<double,double>,std::pair<double,double>> deuterionpair;

    // binned method (worse performance overall)
    std::pair<double,double> pionpair_binned;
    std::pair<double,double> kaonpair_binned;
    std::pair<double,double> protonpair_binned;

    // go over events
    std::cout<<HIFtree_->GetEntries()<<std::endl;
    std::cout<<HIFtreegen_->GetEntries()<<std::endl;

    for(int z=0; z<HIFtree_->GetEntries(); z++){
        HIFtree_->GetEntry(z);
        HIFtreegen_->GetEntry(z);
        
        // fill number of gen particles manually as the field is unknown
        hifgenev_.gen_ntrk = hifgenev_.gen_trk_pt->size();


        ev_.event = z;
        ev_.weight = 1;
        ev_.ntrk = 0;
        
        // loop over pf candidates
        for(int j=0; j<hifev_.ntrk; j++){
          // stop computation if maxtracks limit is hit
          if(ev_.MAXTRACKS==ev_.ntrk){
            std::cout << "ERROR: number of reconstructed tracks reach the maximum of MAXTRACKS =  "<<ev_.MAXTRACKS<<", the ntrk loop is terminated"<<std::endl;
            std::cout <<"\t\t... consider increasing MAXTRACKS !!!"<<std::endl;
            break;
          }

          // don't fill with garbage
          if((hifev_.trk_eta)->at(j)<-2.5 || (hifev_.trk_eta)->at(j)>2.5) continue;

          // fill basic track info
          ev_.trk_p[ev_.ntrk] =  hifev_.trk_pt->at(j)*cosh(hifev_.trk_eta->at(j)) ;
          ev_.trk_pt[ev_.ntrk] = hifev_.trk_pt->at(j) ;
          ev_.trk_eta[ev_.ntrk] = hifev_.trk_eta->at(j) ;
          ev_.trk_phi[ev_.ntrk] = hifev_.trk_phi->at(j) ;
          ev_.trk_q[ev_.ntrk] = hifev_.trk_q->at(j) ;


          ev_.trk_dxy[ev_.ntrk] = hifev_.trk_dxy->at(j) ;
          ev_.trk_dz[ev_.ntrk] = hifev_.trk_dz->at(j) ;
          ev_.trk_numberOfPixelHits[ev_.ntrk] = hifev_.trk_numberOfPixelHits->at(j) ;
          ev_.trk_numberOfHits[ev_.ntrk] = hifev_.trk_numberOfHits->at(j) ;

          ev_.trk_dedx[ev_.ntrk] = hifev_.trk_dedx->at(j) ;


          // fill dE/dx band matching info, returns the tight and loose bands in 1 object
          pionpair = get_minmax(ev_.trk_p[ev_.ntrk],"pion");
          kaonpair = get_minmax(ev_.trk_p[ev_.ntrk],"kaon");
          protonpair = get_minmax(ev_.trk_p[ev_.ntrk],"proton");
          deuterionpair = get_minmax(ev_.trk_p[ev_.ntrk],"deuterion");

          pionpair_binned = get_minmax_binned(ev_.trk_p[ev_.ntrk],pionmap);
          kaonpair_binned = get_minmax_binned(ev_.trk_p[ev_.ntrk],kaonmap);
          protonpair_binned = get_minmax_binned(ev_.trk_p[ev_.ntrk],protonmap);



          //tight bands
          if (hifev_.trk_dedx->at(j) > pionpair.first.first && hifev_.trk_dedx->at(j) < pionpair.first.second){
            ev_.trk_isPi[ev_.ntrk] = 1;
          }
          else{ev_.trk_isPi[ev_.ntrk] = 0;}


          if (hifev_.trk_dedx->at(j) > kaonpair.first.first && hifev_.trk_dedx->at(j) < kaonpair.first.second){
            ev_.trk_isK[ev_.ntrk] = 1;
          }
          else{ev_.trk_isK[ev_.ntrk] = 0;}

          if (hifev_.trk_dedx->at(j) > protonpair.first.first && hifev_.trk_dedx->at(j) < protonpair.first.second){
            ev_.trk_isP[ev_.ntrk] = 1;
          }
          else{ev_.trk_isP[ev_.ntrk] = 0;}

          if (hifev_.trk_dedx->at(j) > deuterionpair.first.first && hifev_.trk_dedx->at(j) < deuterionpair.first.second){
            ev_.trk_isD[ev_.ntrk] = 1;
          }
          else{ev_.trk_isD[ev_.ntrk] = 0;}


          //loose bands
          if (hifev_.trk_dedx->at(j) > pionpair.second.first && hifev_.trk_dedx->at(j) < pionpair.second.second){
            ev_.trk_isPi_loose[ev_.ntrk] = 1;
          }
          else{ev_.trk_isPi_loose[ev_.ntrk] = 0;}


          if (hifev_.trk_dedx->at(j) > kaonpair.second.first && hifev_.trk_dedx->at(j) < kaonpair.second.second){
            ev_.trk_isK_loose[ev_.ntrk] = 1;
          }
          else{ev_.trk_isK_loose[ev_.ntrk] = 0;}

          if (hifev_.trk_dedx->at(j) > protonpair.second.first && hifev_.trk_dedx->at(j) < protonpair.second.second){
            ev_.trk_isP_loose[ev_.ntrk] = 1;
          }
          else{ev_.trk_isP_loose[ev_.ntrk] = 0;}

          if (hifev_.trk_dedx->at(j) > deuterionpair.second.first && hifev_.trk_dedx->at(j) < deuterionpair.second.second){
            ev_.trk_isD_loose[ev_.ntrk] = 1;
          }
          else{ev_.trk_isD_loose[ev_.ntrk] = 0;}


          //similar output for binned method
          if (hifev_.trk_dedx->at(j) > pionpair_binned.first && hifev_.trk_dedx->at(j) < pionpair_binned.second){
            ev_.trk_isPi_binned[ev_.ntrk] = 1;
          }
          else{ev_.trk_isPi_binned[ev_.ntrk] = 0;}

          if (hifev_.trk_dedx->at(j) > kaonpair_binned.first && hifev_.trk_dedx->at(j) < kaonpair_binned.second){
            ev_.trk_isK_binned[ev_.ntrk] = 1;
          }
          else{ev_.trk_isK_binned[ev_.ntrk] = 0;}

          if (hifev_.trk_dedx->at(j) > protonpair_binned.first && hifev_.trk_dedx->at(j) < protonpair_binned.second){
            ev_.trk_isP_binned[ev_.ntrk] = 1;
          }
          else{ev_.trk_isP_binned[ev_.ntrk] = 0;}


          //now look for gen matching information
          dRmatch = 99.;
          pmatch = 99.;
          matchedPdgId=0;

          // loop over gen particles
          for (size_t q = 0; q < hifgenev_.gen_ntrk; ++q){
              // exclude photons, neutrinos, kaons
              if(std::find(v.begin(), v.end(),  hifgenev_.gen_trk_id->at(q)) != v.end() ) continue;
              // if particle is close
              if (deltaR(hifgenev_.gen_trk_eta->at(q), ev_.trk_eta[ev_.ntrk], hifgenev_.gen_trk_phi->at(q), ev_.trk_phi[ev_.ntrk]) < dRmatch){
                  // if dRmatch is already below a certain value, we should decide based on pt agreement as well
                  if (dRmatch<0.1){
                    if ( deltaR(hifgenev_.gen_trk_eta->at(q), ev_.trk_eta[ev_.ntrk], hifgenev_.gen_trk_phi->at(q), ev_.trk_phi[ev_.ntrk]) + 0.5*std::abs(hifgenev_.gen_trk_pt->at(q) - ev_.trk_pt[ev_.ntrk]) < dRmatch + 0.5*std::abs(pmatch - ev_.trk_pt[ev_.ntrk]) ){
                      dRmatch=deltaR(hifgenev_.gen_trk_eta->at(q), ev_.trk_eta[ev_.ntrk], hifgenev_.gen_trk_phi->at(q), ev_.trk_phi[ev_.ntrk]);
                      matchedPdgId=hifgenev_.gen_trk_id->at(q);
                      pmatch=hifgenev_.gen_trk_pt->at(q);
                    }
                  }
                  // else just fill
                  else{
                    dRmatch=deltaR(hifgenev_.gen_trk_eta->at(q), ev_.trk_eta[ev_.ntrk], hifgenev_.gen_trk_phi->at(q), ev_.trk_phi[ev_.ntrk]);
                    matchedPdgId=hifgenev_.gen_trk_id->at(q);
                    pmatch=hifgenev_.gen_trk_pt->at(q);
              }
            }
          }
          // fill gen matching info
          ev_.trk_matchedPdgId[ev_.ntrk] = matchedPdgId;
          ev_.trk_dRmatch[ev_.ntrk] = dRmatch;
          ev_.trk_genPt[ev_.ntrk] = pmatch;

          // only write protons and deuterium 
          ev_.ntrk++;
          
        }
     

    // GEN particles: general info and tracking efficiency info
          ev_.gen_ntrk = 0;
          ev_.gen_nD = 0;
         
          for (int i=0; i<hifgenev_.gen_ntrk; i++){ 
          // check that the total number of reco tracks does not exceed the MAXGENTRACKS
          if(ev_.MAXGENTRACKS==ev_.gen_ntrk){
            std::cout << "ERROR: number of reconstructed tracks reach the maximum of MAXGENTRACKS =  "<<ev_.MAXGENTRACKS<<", the gen_ntrk loop is terminated"<<std::endl;
            std::cout <<"\t\t... consider increasing MAXGENTRACKS !!!"<<std::endl;
            break;
          }

          if(hifgenev_.gen_trk_eta->at(i)<-2.5 || hifgenev_.gen_trk_eta->at(i)>2.5) continue;

          if(hifgenev_.gen_trk_charge->at(i) == 0) continue;
          if(std::find(v.begin(), v.end(),  hifgenev_.gen_trk_id->at(i)) != v.end() ) continue;

          ev_.gen_nD++;
          // write gen info
          ev_.gen_trk_pt[ev_.gen_ntrk] = hifgenev_.gen_trk_pt->at(i);
          ev_.gen_trk_p[ev_.gen_ntrk] = hifgenev_.gen_trk_pt->at(i)*cosh(hifgenev_.gen_trk_eta->at(i));
          ev_.gen_trk_E[ev_.gen_ntrk] = 0.;
          ev_.gen_trk_eta[ev_.gen_ntrk] = hifgenev_.gen_trk_eta->at(i);
          ev_.gen_trk_phi[ev_.gen_ntrk] = hifgenev_.gen_trk_phi->at(i);
          ev_.gen_trk_id[ev_.gen_ntrk] = hifgenev_.gen_trk_id->at(i);

          hasReco = false;
          // check if track matches
          for(int h=0; h<hifev_.ntrk; h++){
                if (deltaR(ev_.gen_trk_eta[ev_.gen_ntrk],hifev_.trk_eta->at(h) ,ev_.gen_trk_phi[ev_.gen_ntrk],hifev_.trk_phi->at(h)) < 0.1){
                        hasReco = true;
                        continue;
                        }
                }
          if (hasReco){ ev_.gen_trk_hasReco[ev_.gen_ntrk] = 1;}
          else{ ev_.gen_trk_hasReco[ev_.gen_ntrk] = 0;}
    
      ev_.gen_ntrk++;

      }   //end of gen loop

      tree_->Fill();
    }        //end of event loop

    //write everything in the output file

    dir->cd();

    tree_->SetName("tree");
    tree_->Write();

    of.Close();
    
    std::cerr << "###done###" << std::endl;
}
