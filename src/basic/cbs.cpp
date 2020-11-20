#include <array>
#include <vector>
#include <algorithm>
#include <math.h>
#include "value.h"
#include "sequence.h"
#include "score_matrix.h"
#include "cbs.h"
#include "config.h"

using std::vector;

#define BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT    (1.e-5)
#define BLAST_KARLIN_LAMBDA_ITER_DEFAULT        17
#define COMPO_SCORE_MIN (-128)
#define LambdaRatioLowerBound 0.5
const int ALPH_TO_NCBI[] = { 1, 16, 13, 4, 3, 15, 5, 7, 8, 9, 11, 10, 12, 6, 14, 17, 18, 20, 22, 19 };

typedef struct Blast_ScoreFreq {
    int         score_min; /**< lowest allowed scores */
    int         score_max; /**< highest allowed scores */
    int         obs_min;   /**< lowest observed (actual) scores */
    int         obs_max;   /**< highest observed (actual) scores */
    double       score_avg; /**< average score, must be negative for local alignment. */
    double* sprob0;    /**< arrays for frequency of given score */
    double* sprob;     /**< arrays for frequency of given score, shifted down by score_min. */
} Blast_ScoreFreq;

/** Underlying frequency ratios for BLOSUM62 */
const double BLOSUM62_FREQRATIOS[28][28] =
{ {0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
 {0.00000000e+00, 3.90294070e+00, 5.64459671e-01, 8.67987664e-01,
  5.44605275e-01, 7.41264113e-01, 4.64893827e-01, 1.05686961e+00,
  5.69364849e-01, 6.32481035e-01, 7.75390239e-01, 6.01945975e-01,
  7.23150342e-01, 5.88307640e-01, 7.54121369e-01, 7.56803943e-01,
  6.12698600e-01, 1.47210399e+00, 9.84401956e-01, 9.36458396e-01,
  4.16548781e-01, 7.50000000e-01, 5.42611869e-01, 7.47274948e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 6.14377313e-01},
 {0.00000000e+00, 5.64459671e-01, 4.43758048e+00, 3.45226274e-01,
  4.74290926e+00, 1.33503378e+00, 3.24101420e-01, 7.38524318e-01,
  9.25449581e-01, 3.33981361e-01, 8.54849426e-01, 2.97257620e-01,
  4.04640322e-01, 4.07083696e+00, 5.53838329e-01, 9.44103648e-01,
  7.02873767e-01, 1.05798620e+00, 8.26250098e-01, 3.51280513e-01,
  2.52855433e-01, 7.50000000e-01, 4.09444638e-01, 1.18382127e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.12208474e-01},
 {0.00000000e+00, 8.67987664e-01, 3.45226274e-01, 1.95765857e+01,
  3.01454345e-01, 2.85934574e-01, 4.38990118e-01, 4.20387870e-01,
  3.55049505e-01, 6.53458801e-01, 3.49128465e-01, 6.42275633e-01,
  6.11354340e-01, 3.97802620e-01, 3.79562691e-01, 3.65781531e-01,
  3.08939296e-01, 7.38415701e-01, 7.40551692e-01, 7.55844055e-01,
  4.49983903e-01, 7.50000000e-01, 4.34203398e-01, 3.16819526e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 6.46828489e-01},
 {0.00000000e+00, 5.44605275e-01, 4.74290926e+00, 3.01454345e-01,
  7.39792738e+00, 1.68781075e+00, 2.98969081e-01, 6.34301019e-01,
  6.78558839e-01, 3.39015407e-01, 7.84090406e-01, 2.86613046e-01,
  3.46454634e-01, 1.55385281e+00, 5.98716826e-01, 8.97081129e-01,
  5.73200024e-01, 9.13504624e-01, 6.94789868e-01, 3.36500142e-01,
  2.32102315e-01, 7.50000000e-01, 3.45683565e-01, 1.38195506e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.07946931e-01},
 {0.00000000e+00, 7.41264113e-01, 1.33503378e+00, 2.85934574e-01,
  1.68781075e+00, 5.46952608e+00, 3.30743991e-01, 4.81267655e-01,
  9.60040718e-01, 3.30522558e-01, 1.30827885e+00, 3.72873704e-01,
  5.00342289e-01, 9.11298183e-01, 6.79202587e-01, 1.90173784e+00,
  9.60797602e-01, 9.50357185e-01, 7.41425610e-01, 4.28943130e-01,
  3.74300212e-01, 7.50000000e-01, 4.96467354e-01, 4.08949895e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.55631838e-01},
 {0.00000000e+00, 4.64893827e-01, 3.24101420e-01, 4.38990118e-01,
  2.98969081e-01, 3.30743991e-01, 8.12879702e+00, 3.40640908e-01,
  6.51990521e-01, 9.45769883e-01, 3.44043119e-01, 1.15459749e+00,
  1.00437163e+00, 3.54288952e-01, 2.87444758e-01, 3.33972402e-01,
  3.80726330e-01, 4.39973597e-01, 4.81693683e-01, 7.45089738e-01,
  1.37437942e+00, 7.50000000e-01, 2.76938063e+00, 3.31992746e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 1.06958025e+00},
 {0.00000000e+00, 1.05686961e+00, 7.38524318e-01, 4.20387870e-01,
  6.34301019e-01, 4.81267655e-01, 3.40640908e-01, 6.87630691e+00,
  4.92966576e-01, 2.75009722e-01, 5.88871736e-01, 2.84504012e-01,
  3.95486600e-01, 8.63711406e-01, 4.77385507e-01, 5.38649627e-01,
  4.49983999e-01, 9.03596525e-01, 5.79271582e-01, 3.36954912e-01,
  4.21690355e-01, 7.50000000e-01, 3.48714366e-01, 5.03463109e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.80638726e-01},
 {0.00000000e+00, 5.69364849e-01, 9.25449581e-01, 3.55049505e-01,
  6.78558839e-01, 9.60040718e-01, 6.51990521e-01, 4.92966576e-01,
  1.35059997e+01, 3.26288125e-01, 7.78887490e-01, 3.80675486e-01,
  5.84132623e-01, 1.22200067e+00, 4.72879831e-01, 1.16798104e+00,
  9.17048021e-01, 7.36731740e-01, 5.57503254e-01, 3.39447442e-01,
  4.44088955e-01, 7.50000000e-01, 1.79790413e+00, 1.04047242e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.58533474e-01},
 {0.00000000e+00, 6.32481035e-01, 3.33981361e-01, 6.53458801e-01,
  3.39015407e-01, 3.30522558e-01, 9.45769883e-01, 2.75009722e-01,
  3.26288125e-01, 3.99792994e+00, 3.96372934e-01, 1.69443475e+00,
  1.47774450e+00, 3.27934752e-01, 3.84662860e-01, 3.82937802e-01,
  3.54751311e-01, 4.43163582e-01, 7.79816110e-01, 2.41751209e+00,
  4.08874390e-01, 7.50000000e-01, 6.30388931e-01, 3.50796872e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.63222650e+00},
 {0.00000000e+00, 7.75390239e-01, 8.54849426e-01, 3.49128465e-01,
  7.84090406e-01, 1.30827885e+00, 3.44043119e-01, 5.88871736e-01,
  7.78887490e-01, 3.96372934e-01, 4.76433717e+00, 4.28270363e-01,
  6.25302816e-01, 9.39841129e-01, 7.03774479e-01, 1.55432308e+00,
  2.07680867e+00, 9.31919141e-01, 7.92905803e-01, 4.56542720e-01,
  3.58930071e-01, 7.50000000e-01, 5.32179333e-01, 1.40344922e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.15284382e-01},
 {0.00000000e+00, 6.01945975e-01, 2.97257620e-01, 6.42275633e-01,
  2.86613046e-01, 3.72873704e-01, 1.15459749e+00, 2.84504012e-01,
  3.80675486e-01, 1.69443475e+00, 4.28270363e-01, 3.79662137e+00,
  1.99429557e+00, 3.10043276e-01, 3.71121724e-01, 4.77325586e-01,
  4.73919278e-01, 4.28893743e-01, 6.60328975e-01, 1.31423573e+00,
  5.68037074e-01, 7.50000000e-01, 6.92059423e-01, 4.13275887e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.94078574e+00},
 {0.00000000e+00, 7.23150342e-01, 4.04640322e-01, 6.11354340e-01,
  3.46454634e-01, 5.00342289e-01, 1.00437163e+00, 3.95486600e-01,
  5.84132623e-01, 1.47774450e+00, 6.25302816e-01, 1.99429557e+00,
  6.48145121e+00, 4.74529655e-01, 4.23898024e-01, 8.64250293e-01,
  6.22623369e-01, 5.98558924e-01, 7.93801616e-01, 1.26893679e+00,
  6.10296214e-01, 7.50000000e-01, 7.08364628e-01, 6.41102583e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 1.78399892e+00},
 {0.00000000e+00, 5.88307640e-01, 4.07083696e+00, 3.97802620e-01,
  1.55385281e+00, 9.11298183e-01, 3.54288952e-01, 8.63711406e-01,
  1.22200067e+00, 3.27934752e-01, 9.39841129e-01, 3.10043276e-01,
  4.74529655e-01, 7.09409488e+00, 4.99932836e-01, 1.00058442e+00,
  8.58630478e-01, 1.23152924e+00, 9.84152635e-01, 3.69033853e-01,
  2.77782896e-01, 7.50000000e-01, 4.86030806e-01, 9.45834265e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.17327197e-01},
 {0.00000000e+00, 7.54121369e-01, 5.53838329e-01, 3.79562691e-01,
  5.98716826e-01, 6.79202587e-01, 2.87444758e-01, 4.77385507e-01,
  4.72879831e-01, 3.84662860e-01, 7.03774479e-01, 3.71121724e-01,
  4.23898024e-01, 4.99932836e-01, 1.28375437e+01, 6.41280589e-01,
  4.81534905e-01, 7.55503259e-01, 6.88897122e-01, 4.43082984e-01,
  2.81833164e-01, 7.50000000e-01, 3.63521119e-01, 6.64534287e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.76634549e-01},
 {0.00000000e+00, 7.56803943e-01, 9.44103648e-01, 3.65781531e-01,
  8.97081129e-01, 1.90173784e+00, 3.33972402e-01, 5.38649627e-01,
  1.16798104e+00, 3.82937802e-01, 1.55432308e+00, 4.77325586e-01,
  8.64250293e-01, 1.00058442e+00, 6.41280589e-01, 6.24442175e+00,
  1.40579606e+00, 9.65555228e-01, 7.91320741e-01, 4.66777931e-01,
  5.09360272e-01, 7.50000000e-01, 6.11094097e-01, 3.58149606e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.38898727e-01},
 {0.00000000e+00, 6.12698600e-01, 7.02873767e-01, 3.08939296e-01,
  5.73200024e-01, 9.60797602e-01, 3.80726330e-01, 4.49983999e-01,
  9.17048021e-01, 3.54751311e-01, 2.07680867e+00, 4.73919278e-01,
  6.22623369e-01, 8.58630478e-01, 4.81534905e-01, 1.40579606e+00,
  6.66557707e+00, 7.67165633e-01, 6.77754679e-01, 4.20072316e-01,
  3.95102106e-01, 7.50000000e-01, 5.55965425e-01, 1.13292384e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.25403989e-01},
 {0.00000000e+00, 1.47210399e+00, 1.05798620e+00, 7.38415701e-01,
  9.13504624e-01, 9.50357185e-01, 4.39973597e-01, 9.03596525e-01,
  7.36731740e-01, 4.43163582e-01, 9.31919141e-01, 4.28893743e-01,
  5.98558924e-01, 1.23152924e+00, 7.55503259e-01, 9.65555228e-01,
  7.67165633e-01, 3.84284741e+00, 1.61392097e+00, 5.65223766e-01,
  3.85303035e-01, 7.50000000e-01, 5.57520051e-01, 9.56235816e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 4.34703235e-01},
 {0.00000000e+00, 9.84401956e-01, 8.26250098e-01, 7.40551692e-01,
  6.94789868e-01, 7.41425610e-01, 4.81693683e-01, 5.79271582e-01,
  5.57503254e-01, 7.79816110e-01, 7.92905803e-01, 6.60328975e-01,
  7.93801616e-01, 9.84152635e-01, 6.88897122e-01, 7.91320741e-01,
  6.77754679e-01, 1.61392097e+00, 4.83210516e+00, 9.80943005e-01,
  4.30934144e-01, 7.50000000e-01, 5.73156574e-01, 7.60725140e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.08974203e-01},
 {0.00000000e+00, 9.36458396e-01, 3.51280513e-01, 7.55844055e-01,
  3.36500142e-01, 4.28943130e-01, 7.45089738e-01, 3.36954912e-01,
  3.39447442e-01, 2.41751209e+00, 4.56542720e-01, 1.31423573e+00,
  1.26893679e+00, 3.69033853e-01, 4.43082984e-01, 4.66777931e-01,
  4.20072316e-01, 5.65223766e-01, 9.80943005e-01, 3.69215640e+00,
  3.74456332e-01, 7.50000000e-01, 6.58038693e-01, 4.43577702e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 1.76339815e+00},
 {0.00000000e+00, 4.16548781e-01, 2.52855433e-01, 4.49983903e-01,
  2.32102315e-01, 3.74300212e-01, 1.37437942e+00, 4.21690355e-01,
  4.44088955e-01, 4.08874390e-01, 3.58930071e-01, 5.68037074e-01,
  6.10296214e-01, 2.77782896e-01, 2.81833164e-01, 5.09360272e-01,
  3.95102106e-01, 3.85303035e-01, 4.30934144e-01, 3.74456332e-01,
  3.81077833e+01, 7.50000000e-01, 2.10980812e+00, 4.26541694e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 5.03239261e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 5.42611869e-01, 4.09444638e-01, 4.34203398e-01,
  3.45683565e-01, 4.96467354e-01, 2.76938063e+00, 3.48714366e-01,
  1.79790413e+00, 6.30388931e-01, 5.32179333e-01, 6.92059423e-01,
  7.08364628e-01, 4.86030806e-01, 3.63521119e-01, 6.11094097e-01,
  5.55965425e-01, 5.57520051e-01, 5.73156574e-01, 6.58038693e-01,
  2.10980812e+00, 7.50000000e-01, 9.83220341e+00, 5.40805192e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 6.66952325e-01},
 {0.00000000e+00, 7.47274948e-01, 1.18382127e+00, 3.16819526e-01,
  1.38195506e+00, 4.08949895e+00, 3.31992746e-01, 5.03463109e-01,
  1.04047242e+00, 3.50796872e-01, 1.40344922e+00, 4.13275887e-01,
  6.41102583e-01, 9.45834265e-01, 6.64534287e-01, 3.58149606e+00,
  1.13292384e+00, 9.56235816e-01, 7.60725140e-01, 4.43577702e-01,
  4.26541694e-01, 7.50000000e-01, 5.40805192e-01, 3.89300249e+00,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 3.87839626e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 2.50000000e-01, 2.50000000e-01, 2.50000000e-01,
  2.50000000e-01, 1.33300000e+00, 2.50000000e-01, 2.50000000e-01},
 {0.00000000e+00, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 7.50000000e-01, 7.50000000e-01, 7.50000000e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 7.50000000e-01},
 {0.00000000e+00, 6.14377313e-01, 3.12208474e-01, 6.46828489e-01,
  3.07946931e-01, 3.55631838e-01, 1.06958025e+00, 2.80638726e-01,
  3.58533474e-01, 2.63222650e+00, 4.15284382e-01, 2.94078574e+00,
  1.78399892e+00, 3.17327197e-01, 3.76634549e-01, 4.38898727e-01,
  4.25403989e-01, 4.34703235e-01, 7.08974203e-01, 1.76339815e+00,
  5.03239261e-01, 7.50000000e-01, 6.66952325e-01, 3.87839626e-01,
  7.50000000e-01, 2.50000000e-01, 7.50000000e-01, 2.81516607e+00} };

const double BLOSUM62_bg[20] =
  { 7.4216205067993410e-02, 5.1614486141284638e-02, 4.4645808512757915e-02,
   5.3626000838554413e-02, 2.4687457167944848e-02, 3.4259650591416023e-02,
   5.4311925684587502e-02, 7.4146941452644999e-02, 2.6212984805266227e-02,
   6.7917367618953756e-02, 9.8907868497150955e-02, 5.8155682303079680e-02,
   2.4990197579643110e-02, 4.7418459742284751e-02, 3.8538003320306206e-02,
   5.7229029476494421e-02, 5.0891364550287033e-02, 1.3029956129972148e-02,
   3.2281512313758580e-02, 7.2919098205619245e-02 };


inline double**
Nlm_DenseMatrixNew(int nrows,
    int ncols)
{
    int i;             /* iteration index */
    double** mat;     /* the new matrix */

    mat = (double**)calloc(nrows, sizeof(double*));
    if (mat != NULL) {
        mat[0] = (double*)malloc((size_t)nrows *
            (size_t)ncols * sizeof(double));
        if (mat[0] != NULL) {
            for (i = 1; i < nrows; i++) {
                mat[i] = &mat[0][i * ncols];
            }
        }
        else {
            free(mat);
            mat = NULL;
        }
    }
    return mat;
}

inline void
Nlm_DenseMatrixFree(double*** mat)
{
    if (*mat != NULL) {
        free((*mat)[0]);
        free(*mat);
    }
    *mat = NULL;
}

inline int BLAST_Gcd(int a, int b)
{
    int   c;

    b = std::abs(b);
    if (b > a)
        c = a, a = b, b = c;

    while (b != 0) {
        c = a % b;
        a = b;
        b = c;
    }
    return a;
}

inline static double
NlmKarlinLambdaNR(double* probs, int d, int low, int high, double lambda0,
    double tolx, int itmax, int maxNewton, int* itn)
{
    int k;
    double x0, x, a = 0, b = 1;
    double f = 4;  /* Larger than any possible value of the poly in [0,1] */
    int isNewton = 0; /* we haven't yet taken a Newton step. */

    assert(d > 0);

    x0 = exp(-lambda0);
    x = (0 < x0 && x0 < 1) ? x0 : .5;

    for (k = 0; k < itmax; k++) { /* all iteration indices k */
        int i;
        double g, fold = f;
        int wasNewton = isNewton; /* If true, then the previous step was a */
                                  /* Newton step */
        isNewton = 0;            /* Assume that this step is not */

        /* Horner's rule for evaluating a polynomial and its derivative */
        g = 0;
        f = probs[low];
        for (i = low + d; i < 0; i += d) {
            g = x * g + f;
            f = f * x + probs[i];
        }
        g = x * g + f;
        f = f * x + probs[0] - 1;
        for (i = d; i <= high; i += d) {
            g = x * g + f;
            f = f * x + probs[i];
        }
        /* End Horner's rule */

        if (f > 0) {
            a = x; /* move the left endpoint */
        }
        else if (f < 0) {
            b = x; /* move the right endpoint */
        }
        else { /* f == 0 */
            break; /* x is an exact solution */
        }
        if (b - a < 2 * a * (1 - b) * tolx) {
            /* The midpoint of the interval converged */
            x = (a + b) / 2; break;
        }

        if (k >= maxNewton ||
            /* If convergence of Newton's method appears to be failing; or */
            (wasNewton && fabs(f) > .9 * fabs(fold)) ||
            /* if the previous iteration was a Newton step but didn't decrease
             * f sufficiently; or */
            g >= 0
            /* if a Newton step will move us away from the desired solution */
            ) { /* then */
          /* bisect */
            x = (a + b) / 2;
        }
        else {
            /* try a Newton step */
            double p = -f / g;
            double y = x + p;
            if (y <= a || y >= b) { /* The proposed iterate is not in (a,b) */
                x = (a + b) / 2;
            }
            else { /* The proposed iterate is in (a,b). Accept it. */
                isNewton = 1;
                x = y;
                if (fabs(p) < tolx * x * (1 - x)) break; /* Converged */
            } /* else the proposed iterate is in (a,b) */
        } /* else try a Newton step. */
    } /* end for all iteration indices k */
    *itn = k;
    return -log(x) / d;
}

inline double
Blast_KarlinLambdaNR(Blast_ScoreFreq* sfp, double initialLambdaGuess)
{
    int  low;        /* Lowest score (must be negative)  */
    int  high;       /* Highest score (must be positive) */
    int     itn;
    int i, d;
    double* sprob;
    double   returnValue;

    low = sfp->obs_min;
    high = sfp->obs_max;
    if (sfp->score_avg >= 0.) {   /* Expected score must be negative */
        return -1.0;
    }
    //if (BlastScoreChk(low, high) != 0) return -1.;

    sprob = sfp->sprob;
    /* Find greatest common divisor of all scores */
    for (i = 1, d = -low; i <= high - low && d > 1; ++i) {
        if (sprob[i + low] != 0.0) {
            d = BLAST_Gcd(d, i);
        }
    }
    returnValue =
        NlmKarlinLambdaNR(sprob, d, low, high,
            initialLambdaGuess,
            BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT,
            20, 20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT, &itn);


    return returnValue;
}

inline static double
s_CalcLambda(double probs[], int min_score, int max_score, double lambda0)
{

    int i;                 /* loop index */
    int score_range;       /* range of possible scores */
    double avg;            /* expected score of aligning two characters */
    Blast_ScoreFreq freq;  /* score frequency data */

    score_range = max_score - min_score + 1;
    avg = 0.0;
    for (i = 0; i < score_range; i++) {
        avg += (min_score + i) * probs[i];
    }
    freq.score_min = min_score;
    freq.score_max = max_score;
    freq.obs_min = min_score;
    freq.obs_max = max_score;
    freq.sprob0 = probs;
    freq.sprob = &probs[-min_score];
    freq.score_avg = avg;

    return Blast_KarlinLambdaNR(&freq, lambda0);
}

inline static void s_GetScoreRange(int* obs_min, int* obs_max,
    const int* const* matrix, int rows)
{
    int aa;                    /* index of an amino-acid in the 20
                                  letter alphabet */
    int irow, jcol;            /* matrix row and column indices */
    int minScore, maxScore;    /* largest and smallest observed scores */

    minScore = maxScore = 0;
    for (irow = 0; irow < rows; irow++) {
        for (aa = 0; aa < 20; aa++) {
            jcol = aa;
            if (matrix[irow][jcol] < minScore)
                minScore = matrix[irow][jcol];
            if (matrix[irow][jcol] > maxScore)
                maxScore = matrix[irow][jcol];
        }
    }
    *obs_min = minScore;
    *obs_max = maxScore;
}

inline static int
s_GetMatrixScoreProbs(double** scoreProb, int* obs_min, int* obs_max,
    const int* const* matrix, int alphsize,
    const double* subjectProbArray,
    const double* queryProbArray)
{
    int aa;          /* index of an amino-acid in the 20 letter
                        alphabet */
    int irow, jcol;  /* matrix row and column indices */
    double* sprob;  /* a pointer to the element of the score
                        probabilities array that represents the
                        probability of the score 0*/
    int minScore;    /* smallest score in matrix; the same value as
                        (*obs_min). */
    int range;       /* the range of scores in the matrix */

    s_GetScoreRange(obs_min, obs_max, matrix, alphsize);
    minScore = *obs_min;
    range = *obs_max - *obs_min + 1;
    *scoreProb = (double*)calloc(range, sizeof(double));
    if (*scoreProb == NULL) {
        return -1;
    }
    sprob = &((*scoreProb)[-(*obs_min)]); /*center around 0*/
    for (irow = 0; irow < alphsize; irow++) {
        for (aa = 0; aa < 20; aa++) {
            jcol = aa;
            if (matrix[irow][jcol] >= minScore) {
                sprob[matrix[irow][jcol]] +=
                    (queryProbArray[irow] * subjectProbArray[jcol]);
            }
        }
    }
    return 0;
}

inline void
Blast_FreqRatioToScore(double** matrix, int rows, int cols, double Lambda)
{
    int i;
    for (i = 0; i < rows; i++) {
        int j;
        for (j = 0; j < cols; j++) {
            if (0.0 == matrix[i][j]) {
                matrix[i][j] = COMPO_SCORE_MIN;
            }
            else {
                matrix[i][j] = log(matrix[i][j]) / Lambda;
            }
        }
    }
}

inline static void
s_RoundScoreMatrix(int** matrix, int rows, int cols,
    double** floatScoreMatrix)
{
    int p, c; /*indices over positions and characters*/

    for (p = 0; p < rows; p++) {
        for (c = 0; c < cols; c++) {
            if (floatScoreMatrix[p][c] < INT_MIN) {
                matrix[p][c] = INT_MIN;
            }
            else {
                matrix[p][c] = std::round(floatScoreMatrix[p][c]);
            }
        }
    }
}

inline static int
s_ScaleSquareMatrix(int** matrix, int alphsize,
    const double row_prob[], const double col_prob[],
    double Lambda)
{
    double** scores;     /* a double precision matrix of scores */
    int i;                /* iteration index */

    scores = Nlm_DenseMatrixNew(alphsize, alphsize);
    if (scores == 0) return -1;

    for (i = 0; i < alphsize; i++) {
        for (int j = 0; j < alphsize; ++j)
            scores[i][j] = BLOSUM62_FREQRATIOS[ALPH_TO_NCBI[i]][ALPH_TO_NCBI[j]];
        //memcpy(scores[i], freq_ratios[i], alphsize * sizeof(double));
    }
    Blast_FreqRatioToScore(scores, alphsize, alphsize, Lambda);
    //s_SetXUOScores(scores, alphsize, row_prob, col_prob);
    s_RoundScoreMatrix(matrix, alphsize, alphsize, scores);
    /*for (i = 0; i < alphsize; i++) {
        matrix[i][(int)STOP_LETTER] = start_matrix[i][(int)STOP_LETTER];
        matrix[(int)STOP_LETTER][i] = start_matrix[(int)STOP_LETTER][i];
    }*/
    Nlm_DenseMatrixFree(&scores);

    return 0;
}

inline int
Blast_CompositionBasedStats(int** matrix, double* LambdaRatio,
    const int* const* matrix_in,
    const double queryProb[], const double resProb[])
{
    double correctUngappedLambda; /* new value of ungapped lambda */
    int obs_min, obs_max;         /* smallest and largest score in the
                                     unscaled matrix */
    double* scoreArray;           /* an array of score probabilities */
    int out_of_memory;            /* status flag to indicate out of memory */

    /*if (ungappedLambda == 0.0) {

        s_GetMatrixScoreProbs(&scoreArray, &obs_min, &obs_max, matrix_in, 20, BLOSUM62_bg, BLOSUM62_bg);
        ungappedLambda = s_CalcLambda(scoreArray, obs_min, obs_max, 0.5);
        //std::cerr << "lambda=" << ungappedLambda << endl;
        free(scoreArray);

    }*/

    out_of_memory = s_GetMatrixScoreProbs(&scoreArray, &obs_min, &obs_max, matrix_in, 20, resProb, queryProb);
    const double ungappedLambda = BLOSUM62_UNGAPPED_LAMBDA / config.cbs_matrix_scale;

    if (out_of_memory)
        return -1;
    correctUngappedLambda =
        s_CalcLambda(scoreArray, obs_min, obs_max, ungappedLambda);

    /* calc_lambda will return -1 in the case where the
     * expected score is >=0; however, because of the MAX statement 3
     * lines below, LambdaRatio should always be > 0; the succeeding
     * test is retained as a vestige, in case one wishes to remove the
     * MAX statement and allow LambdaRatio to take on the error value
     * -1 */
    *LambdaRatio = correctUngappedLambda / ungappedLambda;
    //if (0 == pValueAdjustment)
    *LambdaRatio = std::min(1.0, *LambdaRatio);
    *LambdaRatio = std::max(*LambdaRatio, LambdaRatioLowerBound);

    if (*LambdaRatio > 0) {
        double scaledLambda = ungappedLambda / (*LambdaRatio);

        s_ScaleSquareMatrix(matrix, 20,
            queryProb, resProb, scaledLambda);

    }
    free(scoreArray);

    return 0;
}

std::array<double, 20> composition(const sequence& s) {
    std::array<double, 20> r;
    r.fill(0.0);
    int n = 0;
    for (size_t i = 0; i < s.length(); ++i) {
        int l = s[i];
        if (l < 20) {
            ++r[l];
            ++n;
        }
    }
    for (int i = 0; i < 20; ++i)
        r[i] /= n;
    return r;
}

TargetMatrix::TargetMatrix(const double* query_comp, const sequence& target) :
    scores(32 * AMINO_ACID_COUNT),
    scores32(32 * AMINO_ACID_COUNT)
{
    const vector<const int*> p1 = score_matrix.matrix32_scaled_pointers();

    auto c = composition(target);
    std::vector<int> s2(20 * 20);
    std::vector<int*> p2(20);
    for (int i = 0; i < 20; ++i) {
        p2[i] = &s2[i * 20];
    }
    Blast_CompositionBasedStats(p2.data(), &lambda_ratio, p1.data(), query_comp, c.data());
    for (size_t i = 0; i < AMINO_ACID_COUNT; ++i) {
        for (size_t j = 0; j < AMINO_ACID_COUNT; ++j)
            if (i < 20 && j < 20) {
                scores[i * 32 + j] = s2[i * 20 + j];
                scores32[i * 32 + j] = s2[i * 20 + j];
                //std::cerr << s2[i * 20 + j] << ' ';
            }
            else {
                scores[i * 32 + j] = score_matrix(i, j) * config.cbs_matrix_scale;
                scores32[i * 32 + j] = score_matrix(i, j) * config.cbs_matrix_scale;
            }
        //std::cerr << endl;
    }
}