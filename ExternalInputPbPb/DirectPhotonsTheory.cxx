/*
 *  DirectPhotonsTheory.cxx
 *  
 *  Created by Sami Räsänen on 9th February 2014
 *  Contact: sami.s.rasanen@jyu.fi
 *
 */
#include <DirectPhotonsTheory.h>

DirectPhotonsTheory::DirectPhotonsTheory(void){
	
    // =========================
    // The results from Holopainen, Räsänen, Eskola; Phys.Rev. C84 (2011) 064903
    // =========================
	const int nPointsHolopainen = 56;
    
	const double pTHolopainen[nPointsHolopainen] = {
        0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4,
        1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4,
        2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4,
        3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4,
        4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4,
        5.5, 5.6, 5.7, 5.8, 5.9, 6.0 };
	
	// ===================
	// The default results - lattice Equation of State (EoS),
    // swithching temperature of rates Tc = 170 MeV,
    // using Turbide-Rapp-Gale (TRG) hadron gas rates
	// ===================

    // 0-20 % centrality, default
	const double thermalHolopainen00to20[nPointsHolopainen] = {
        3.9346e+00, 2.5044e+00, 1.6633e+00, 1.1377e+00, 7.9517e-01, 5.6509e-01,
        4.0700e-01, 2.9644e-01, 2.1802e-01, 1.6173e-01, 1.2092e-01, 9.1070e-02,
        6.9055e-02, 5.2700e-02, 4.0465e-02, 3.1253e-02, 2.4273e-02, 1.8953e-02,
        1.4875e-02, 1.1731e-02, 9.2954e-03, 7.3981e-03, 5.9130e-03, 4.7451e-03,
        3.8223e-03, 3.0901e-03, 2.5067e-03, 2.0400e-03, 1.6652e-03, 1.3631e-03,
        1.1189e-03, 9.2072e-04, 7.5947e-04, 6.2788e-04, 5.2019e-04, 4.3183e-04,
        3.5916e-04, 2.9925e-04, 2.4976e-04, 2.0878e-04, 1.7479e-04, 1.4655e-04,
        1.2304e-04, 1.0343e-04, 8.7055e-05, 7.3359e-05, 6.1886e-05, 5.2264e-05,
        4.4183e-05, 3.7388e-05, 3.1668e-05, 2.6847e-05, 2.2779e-05, 1.9344e-05,
        1.6440e-05, 1.3983e-05 };
    
    // 0-20 % centrality, default
	const double v2Holopainen00to20[nPointsHolopainen] = {
        9.9374e-03, 1.1053e-02, 1.2383e-02, 1.3802e-02, 1.5214e-02, 1.6553e-02,
        1.7777e-02, 1.8861e-02, 1.9790e-02, 2.0557e-02, 2.1159e-02, 2.1597e-02,
        2.1874e-02, 2.1998e-02, 2.1974e-02, 2.1814e-02, 2.1528e-02, 2.1129e-02,
        2.0630e-02, 2.0044e-02, 1.9386e-02, 1.8668e-02, 1.7905e-02, 1.7107e-02,
        1.6288e-02, 1.5456e-02, 1.4622e-02, 1.3793e-02, 1.2976e-02, 1.2178e-02,
        1.1402e-02, 1.0653e-02, 9.9334e-03, 9.2455e-03, 8.5906e-03, 7.9695e-03,
        7.3825e-03, 6.8295e-03, 6.3101e-03, 5.8235e-03, 5.3688e-03, 4.9447e-03,
        4.5501e-03, 4.1835e-03, 3.8436e-03, 3.5289e-03, 3.2380e-03, 2.9694e-03,
        2.7217e-03, 2.4937e-03, 2.2838e-03, 2.0910e-03, 1.9139e-03, 1.7515e-03,
        1.6026e-03, 1.4662e-03 };

    // 20-40 % centrality, default
	const double thermalHolopainen20to40[nPointsHolopainen] = {
        1.3419e+00, 8.5228e-01, 5.6493e-01, 3.8567e-01, 2.6901e-01, 1.9073e-01,
        1.3700e-01, 9.9476e-02, 7.2900e-02, 5.3861e-02, 4.0088e-02, 3.0039e-02,
        2.2651e-02, 1.7183e-02, 1.3108e-02, 1.0054e-02, 7.7516e-03, 6.0061e-03,
        4.6759e-03, 3.6570e-03, 2.8727e-03, 2.2661e-03, 1.7947e-03, 1.4269e-03,
        1.1385e-03, 9.1161e-04, 7.3232e-04, 5.9012e-04, 4.7693e-04, 3.8653e-04,
        3.1409e-04, 2.5586e-04, 2.0892e-04, 1.7097e-04, 1.4020e-04, 1.1520e-04,
        9.4831e-05, 7.8203e-05, 6.4598e-05, 5.3445e-05, 4.4283e-05, 3.6744e-05,
        3.0529e-05, 2.5398e-05, 2.1155e-05, 1.7642e-05, 1.4728e-05, 1.2308e-05,
        1.0296e-05, 8.6215e-06, 7.2258e-06, 6.0614e-06, 5.0889e-06, 4.2760e-06,
        3.5957e-06, 3.0258e-06 };

    // 20-40 % centrality, default
	const double v2Holopainen20to40[nPointsHolopainen] = {
        2.8489e-02, 3.1835e-02, 3.5822e-02, 4.0067e-02, 4.4274e-02, 4.8246e-02,
        5.1861e-02, 5.5049e-02, 5.7776e-02, 6.0025e-02, 6.1795e-02, 6.3092e-02,
        6.3931e-02, 6.4330e-02, 6.4313e-02, 6.3907e-02, 6.3144e-02, 6.2058e-02,
        6.0684e-02, 5.9059e-02, 5.7222e-02, 5.5208e-02, 5.3055e-02, 5.0796e-02,
        4.8464e-02, 4.6088e-02, 4.3695e-02, 4.1308e-02, 3.8947e-02, 3.6631e-02,
        3.4373e-02, 3.2186e-02, 3.0078e-02, 2.8056e-02, 2.6126e-02, 2.4291e-02,
        2.2552e-02, 2.0910e-02, 1.9363e-02, 1.7911e-02, 1.6551e-02, 1.5280e-02,
        1.4095e-02, 1.2991e-02, 1.1966e-02, 1.1015e-02, 1.0135e-02, 9.3198e-03,
        8.5671e-03, 7.8725e-03, 7.2323e-03, 6.6428e-03, 6.1005e-03, 5.6019e-03,
        5.1439e-03, 4.7236e-03 };
    
    // 0-40 % centrality, default
	const double thermalHolopainen00to40[nPointsHolopainen] = {
		2.5969e+00, 1.6531e+00, 1.0983e+00, 7.5162e-01, 5.2562e-01, 3.7371e-01,
		2.6926e-01, 1.9616e-01, 1.4427e-01, 1.0700e-01, 7.9972e-02, 6.0192e-02,
		4.5603e-02, 3.4765e-02, 2.6660e-02, 2.0560e-02, 1.5941e-02, 1.2424e-02,
		9.7310e-03, 7.6578e-03, 6.0536e-03, 4.8062e-03, 3.8315e-03, 3.0665e-03, 
		2.4634e-03, 1.9858e-03, 1.6062e-03, 1.3033e-03, 1.0607e-03, 8.6565e-04,
		7.0837e-04, 5.8111e-04, 4.7785e-04, 3.9382e-04, 3.2525e-04, 2.6916e-04,
		2.2315e-04, 1.8535e-04, 1.5420e-04, 1.2849e-04, 1.0724e-04, 8.9622e-05,
		7.5003e-05, 6.2851e-05, 5.2732e-05, 4.4294e-05, 3.7248e-05, 3.1356e-05,
		2.6423e-05, 2.2288e-05, 1.8818e-05, 1.5902e-05, 1.3449e-05, 1.1385e-05,
		9.6445e-06, 8.1765e-06 };
	 
    // 0-40 % centrality, default
	const double v2Holopainen00to40[nPointsHolopainen] = {
		1.8535e-02, 2.0629e-02, 2.3119e-02, 2.5770e-02, 2.8403e-02, 3.0896e-02,
		3.3172e-02, 3.5187e-02, 3.6917e-02, 3.8349e-02, 3.9481e-02, 4.0315e-02,
		4.0857e-02, 4.1119e-02, 4.1116e-02, 4.0863e-02, 4.0381e-02, 3.9690e-02,
		3.8815e-02, 3.7779e-02, 3.6606e-02, 3.5320e-02, 3.3944e-02, 3.2501e-02,
		3.1011e-02, 2.9492e-02, 2.7962e-02, 2.6436e-02, 2.4927e-02, 2.3446e-02,
		2.2002e-02, 2.0603e-02, 1.9254e-02, 1.7961e-02, 1.6725e-02, 1.5549e-02,
		1.4435e-02, 1.3382e-02, 1.2391e-02, 1.1459e-02, 1.0586e-02, 9.7698e-03,
		9.0083e-03, 8.2992e-03, 7.6400e-03, 7.0283e-03, 6.4616e-03, 5.9371e-03,
		5.4525e-03, 5.0052e-03, 4.5929e-03, 4.2131e-03, 3.8637e-03, 3.5426e-03,
		3.2476e-03, 2.9769e-03 };
	
	// ===================
    // Lattice EoS and TRG rates, but swithching temperature 200 MeV
	// ===================
    
    // 0-20 % centrality, Tswitch = 200 MeV
	const double thermalHolopainen00to20Tc200[nPointsHolopainen] = {
        2.7314e+00, 1.7692e+00, 1.1985e+00, 8.3659e-01, 5.9647e-01, 4.3201e-01,
        3.1678e-01, 2.3464e-01, 1.7531e-01, 1.3199e-01, 1.0006e-01, 7.6349e-02,
        5.8609e-02, 4.5250e-02, 3.5128e-02, 2.7413e-02, 2.1500e-02, 1.6944e-02,
        1.3414e-02, 1.0667e-02, 8.5171e-03, 6.8278e-03, 5.4941e-03, 4.4367e-03,
        3.5948e-03, 2.9220e-03, 2.3822e-03, 1.9477e-03, 1.5966e-03, 1.3121e-03,
        1.0809e-03, 8.9239e-04, 7.3831e-04, 6.1205e-04, 5.0834e-04, 4.2295e-04,
        3.5250e-04, 2.9425e-04, 2.4600e-04, 2.0596e-04, 1.7267e-04, 1.4495e-04,
        1.2183e-04, 1.0252e-04, 8.6367e-05, 7.2839e-05, 6.1494e-05, 5.1967e-05,
        4.3958e-05, 3.7218e-05, 3.1539e-05, 2.6749e-05, 2.2705e-05, 1.9288e-05,
        1.6398e-05, 1.3950e-05 };
    
    // 0-20 % centrality, Tswitch = 200 MeV
    const double v2Holopainen00to20Tc200[nPointsHolopainen] = {
        6.9200e-03, 6.8717e-03, 7.1218e-03, 7.5906e-03, 8.1902e-03, 8.8486e-03,
        9.5125e-03, 1.0144e-02, 1.0719e-02, 1.1218e-02, 1.1631e-02, 1.1952e-02,
        1.2178e-02, 1.2311e-02, 1.2353e-02, 1.2310e-02, 1.2187e-02, 1.1993e-02,
        1.1736e-02, 1.1424e-02, 1.1067e-02, 1.0673e-02, 1.0250e-02, 9.8059e-03,
        9.3477e-03, 8.8817e-03, 8.4133e-03, 7.9474e-03, 7.4879e-03, 7.0383e-03,
        6.6013e-03, 6.1790e-03, 5.7731e-03, 5.3847e-03, 5.0146e-03, 4.6633e-03,
        4.3309e-03, 4.0174e-03, 3.7226e-03, 3.4460e-03, 3.1870e-03, 2.9452e-03,
        2.7197e-03, 2.5098e-03, 2.3148e-03, 2.1338e-03, 1.9661e-03, 1.8110e-03,
        1.6676e-03, 1.5351e-03, 1.4129e-03, 1.3003e-03, 1.1966e-03, 1.1012e-03,
        1.0135e-03, 9.3280e-04 };

    // 20-40 % centrality, Tswitch = 200 MeV
	const double thermalHolopainen20to40Tc200[nPointsHolopainen] = {
        9.6803e-01, 6.2079e-01, 4.1730e-01, 2.8942e-01, 2.0514e-01, 1.4775e-01,
        1.0774e-01, 7.9355e-02, 5.8945e-02, 4.4111e-02, 3.3231e-02, 2.5189e-02,
        1.9203e-02, 1.4720e-02, 1.1342e-02, 8.7824e-03, 6.8330e-03, 5.3406e-03,
        4.1924e-03, 3.3049e-03, 2.6156e-03, 2.0780e-03, 1.6569e-03, 1.3256e-03,
        1.0641e-03, 8.5674e-04, 6.9183e-04, 5.6020e-04, 4.5480e-04, 3.7013e-04,
        3.0193e-04, 2.4683e-04, 2.0221e-04, 1.6597e-04, 1.3648e-04, 1.1242e-04,
        9.2762e-05, 7.6658e-05, 6.3444e-05, 5.2581e-05, 4.3637e-05, 3.6260e-05,
        3.0167e-05, 2.5126e-05, 2.0951e-05, 1.7488e-05, 1.4613e-05, 1.2221e-05,
        1.0231e-05, 8.5725e-06, 7.1890e-06, 6.0337e-06, 5.0681e-06, 4.2602e-06,
        3.5838e-06, 3.0169e-06 };
    
    // 20-40 % centrality, Tswitch = 200 MeV
    const double v2Holopainen20to40Tc200[nPointsHolopainen] = {
        2.1970e-02, 2.1735e-02, 2.2551e-02, 2.4126e-02, 2.6134e-02, 2.8315e-02,
        3.0480e-02, 3.2508e-02, 3.4319e-02, 3.5864e-02, 3.7116e-02, 3.8065e-02,
        3.8709e-02, 3.9055e-02, 3.9120e-02, 3.8920e-02, 3.8480e-02, 3.7824e-02,
        3.6979e-02, 3.5973e-02, 3.4831e-02, 3.3581e-02, 3.2246e-02, 3.0851e-02,
        2.9416e-02, 2.7959e-02, 2.6498e-02, 2.5046e-02, 2.3616e-02, 2.2217e-02,
        2.0859e-02, 1.9546e-02, 1.8284e-02, 1.7077e-02, 1.5926e-02, 1.4834e-02,
        1.3800e-02, 1.2824e-02, 1.1906e-02, 1.1044e-02, 1.0237e-02, 9.4817e-03,
        8.7772e-03, 8.1207e-03, 7.5101e-03, 6.9428e-03, 6.4165e-03, 5.9287e-03,
        5.4771e-03, 5.0595e-03, 4.6735e-03, 4.3172e-03, 3.9883e-03, 3.6851e-03,
        3.4057e-03, 3.1482e-03 };
    
    // 0-40 % centrality, Tswitch = 200 MeV
	const double thermalHolopainen00to40Tc200[nPointsHolopainen] = {
		1.8977e+00, 1.2200e+00, 8.2193e-01, 5.7129e-01, 4.0585e-01, 2.9302e-01,
		2.1423e-01, 1.5825e-01, 1.1792e-01, 8.8548e-02, 6.6958e-02, 5.0959e-02,
		3.9017e-02, 3.0045e-02, 2.3262e-02, 1.8105e-02, 1.4161e-02, 1.1129e-02,
		8.7856e-03, 6.9661e-03, 5.5463e-03, 4.4332e-03, 3.5568e-03, 2.8637e-03,
		2.3134e-03, 1.8748e-03, 1.5239e-03, 1.2421e-03, 1.0151e-03, 8.3172e-04,
		6.8306e-04, 5.6221e-04, 4.6372e-04, 3.8325e-04, 3.1733e-04, 2.6321e-04,
		2.1870e-04, 1.8200e-04, 1.5168e-04, 1.2660e-04, 1.0581e-04, 8.8547e-05,
		7.4193e-05, 6.2239e-05, 5.2270e-05, 4.3945e-05, 3.6984e-05, 3.1157e-05,
		2.6272e-05, 2.2174e-05, 1.8731e-05, 1.5836e-05, 1.3400e-05, 1.1347e-05,
		9.6158e-06, 8.1548e-06 };
	
    // 0-40 % centrality, Tswitch = 200 MeV
	const double v2Holopainen00to40Tc200[nPointsHolopainen] = {
		1.4195e-02, 1.4011e-02, 1.4473e-02, 1.5407e-02, 1.6615e-02, 1.7938e-02,
		1.9262e-02, 2.0508e-02, 2.1627e-02, 2.2587e-02, 2.3369e-02, 2.3963e-02,
		2.4369e-02, 2.4589e-02, 2.4633e-02, 2.4510e-02, 2.4236e-02, 2.3825e-02,
		2.3294e-02, 2.2661e-02, 2.1942e-02, 2.1154e-02, 2.0312e-02, 1.9431e-02,
		1.8525e-02, 1.7605e-02, 1.6682e-02, 1.5764e-02, 1.4860e-02, 1.3976e-02,
		1.3117e-02, 1.2287e-02, 1.1488e-02, 1.0725e-02, 9.9964e-03, 9.3049e-03,
		8.6504e-03, 8.0327e-03, 7.4514e-03, 6.9055e-03, 6.3942e-03, 5.9161e-03,
		5.4700e-03, 5.0544e-03, 4.6678e-03, 4.3088e-03, 3.9758e-03, 3.6673e-03,
		3.3818e-03, 3.1178e-03, 2.8740e-03, 2.6491e-03, 2.4417e-03, 2.2506e-03,
		2.0746e-03, 1.9127e-03 };
	
	// ==================
    // Bag EoS with critical temperature Tc = 165 MeV and
    // "old" emission rates by Kapusta, et all, at 1992
	// ===================

    // 0-20 % centrality, bag EoS + R92
	const double thermalHolopainen00to20Bag[nPointsHolopainen] = {
        3.2928e+00, 2.1146e+00, 1.4195e+00, 9.8176e-01, 6.9344e-01, 4.9746e-01,
        3.6115e-01, 2.6471e-01, 1.9560e-01, 1.4554e-01, 1.0898e-01, 8.2066e-02,
        6.2133e-02, 4.7282e-02, 3.6156e-02, 2.7779e-02, 2.1441e-02, 1.6623e-02,
        1.2944e-02, 1.0123e-02, 7.9497e-03, 6.2689e-03, 4.9633e-03, 3.9451e-03,
        3.1477e-03, 2.5208e-03, 2.0260e-03, 1.6339e-03, 1.3221e-03, 1.0732e-03,
        8.7383e-04, 7.1355e-04, 5.8428e-04, 4.7969e-04, 3.9480e-04, 3.2570e-04,
        2.6929e-04, 2.2313e-04, 1.8524e-04, 1.5408e-04, 1.2839e-04, 1.0716e-04,
        8.9580e-05, 7.4997e-05, 6.2877e-05, 5.2786e-05, 4.4370e-05, 3.7340e-05,
        3.1459e-05, 2.6533e-05, 2.2402e-05, 1.8931e-05, 1.6013e-05, 1.3557e-05,
        1.1487e-05, 9.7411e-06 };
    
    // 0-20 % centrality, bag EoS + R92
	const double v2Holopainen00to20Bag[nPointsHolopainen] = {
        1.2234e-02, 1.3241e-02, 1.5023e-02, 1.7280e-02, 1.9778e-02, 2.2347e-02,
        2.4875e-02, 2.7286e-02, 2.9531e-02, 3.1576e-02, 3.3396e-02, 3.4974e-02,
        3.6297e-02, 3.7358e-02, 3.8153e-02, 3.8682e-02, 3.8950e-02, 3.8966e-02,
        3.8741e-02, 3.8291e-02, 3.7632e-02, 3.6785e-02, 3.5772e-02, 3.4616e-02,
        3.3340e-02, 3.1967e-02, 3.0520e-02, 2.9020e-02, 2.7488e-02, 2.5943e-02,
        2.4400e-02, 2.2875e-02, 2.1381e-02, 1.9927e-02, 1.8522e-02, 1.7173e-02,
        1.5886e-02, 1.4663e-02, 1.3506e-02, 1.2418e-02, 1.1397e-02, 1.0443e-02,
        9.5550e-03, 8.7301e-03, 7.9662e-03, 7.2605e-03, 6.6102e-03, 6.0121e-03,
        5.4631e-03, 4.9601e-03, 4.5001e-03, 4.0799e-03, 3.6967e-03, 3.3477e-03,
        3.0302e-03, 2.7417e-03 };
    
    // 20-40 % centrality, bag EoS + R92
	const double thermalHolopainen20to40Bag[nPointsHolopainen] = {
        1.1395e+00, 7.3034e-01, 4.8967e-01, 3.3837e-01, 2.3879e-01, 1.7111e-01,
        1.2401e-01, 9.0685e-02, 6.6799e-02, 4.9508e-02, 3.6893e-02, 2.7627e-02,
        2.0784e-02, 1.5703e-02, 1.1914e-02, 9.0755e-03, 6.9410e-03, 5.3294e-03,
        4.1080e-03, 3.1788e-03, 2.4692e-03, 1.9254e-03, 1.5070e-03, 1.1839e-03,
        9.3346e-04, 7.3866e-04, 5.8657e-04, 4.6738e-04, 3.7364e-04, 2.9966e-04,
        2.4106e-04, 1.9449e-04, 1.5736e-04, 1.2767e-04, 1.0384e-04, 8.4659e-05,
        6.9182e-05, 5.6657e-05, 4.6495e-05, 3.8229e-05, 3.1490e-05, 2.5984e-05,
        2.1475e-05, 1.7776e-05, 1.4735e-05, 1.2231e-05, 1.0165e-05, 8.4587e-06,
        7.0467e-06, 5.8767e-06, 4.9060e-06, 4.0996e-06, 3.4289e-06, 2.8705e-06,
        2.4050e-06, 2.0165e-06 };
    
    // 20-40 % centrality, bag EoS + R92
	const double v2Holopainen20to40Bag[nPointsHolopainen] = {
        3.4238e-02, 3.7335e-02, 4.2675e-02, 4.9360e-02, 5.6680e-02, 6.4137e-02,
        7.1407e-02, 7.8285e-02, 8.4644e-02, 9.0402e-02, 9.5504e-02, 9.9915e-02,
        1.0361e-01, 1.0658e-01, 1.0881e-01, 1.1032e-01, 1.1111e-01, 1.1120e-01,
        1.1063e-01, 1.0944e-01, 1.0766e-01, 1.0535e-01, 1.0258e-01, 9.9383e-02,
        9.5839e-02, 9.2008e-02, 8.7951e-02, 8.3730e-02, 7.9401e-02, 7.5018e-02,
        7.0629e-02, 6.6278e-02, 6.2002e-02, 5.7833e-02, 5.3797e-02, 4.9915e-02,
        4.6204e-02, 4.2674e-02, 3.9334e-02, 3.6187e-02, 3.3233e-02, 3.0472e-02,
        2.7900e-02, 2.5510e-02, 2.3297e-02, 2.1252e-02, 1.9367e-02, 1.7633e-02,
        1.6042e-02, 1.4584e-02, 1.3250e-02, 1.2032e-02, 1.0921e-02, 9.9085e-03,
        8.9872e-03, 8.1496e-03 };
    
    // 0-40 % centrality, bag EoS + R92
	const double thermalHolopainen00to40Bag[nPointsHolopainen] = {
		2.1839e+00, 1.4022e+00, 9.4173e-01, 6.5189e-01, 4.6096e-01, 3.3106e-01,
		2.4060e-01, 1.7651e-01, 1.3051e-01, 9.7145e-02, 7.2742e-02, 5.4765e-02,
		4.1438e-02, 3.1504e-02, 2.4060e-02, 1.8456e-02, 1.4218e-02, 1.0999e-02,
		8.5437e-03, 6.6634e-03, 5.2177e-03, 4.1017e-03, 3.2368e-03, 2.5640e-03,
		2.0385e-03, 1.6265e-03, 1.3024e-03, 1.0464e-03, 8.4342e-04, 6.8199e-04,
		5.5313e-04, 4.4991e-04, 3.6698e-04, 3.0012e-04, 2.4605e-04, 2.0221e-04,
		1.6655e-04, 1.3748e-04, 1.1371e-04, 9.4226e-05, 7.8224e-05, 6.5049e-05,
		5.4181e-05, 4.5196e-05, 3.7756e-05, 3.1583e-05, 2.6453e-05, 2.2183e-05,
		1.8624e-05, 1.5652e-05, 1.3169e-05, 1.1090e-05, 9.3481e-06, 7.8867e-06,
		6.6594e-06, 5.6277e-06 };

    // 0-40 % centrality, bag EoS + R92
	const double v2Holopainen00to40Bag[nPointsHolopainen] = {
		2.2763e-02, 2.4652e-02, 2.7968e-02, 3.2153e-02, 3.6769e-02, 4.1506e-02,
		4.6158e-02, 5.0591e-02, 5.4718e-02, 5.8482e-02, 6.1842e-02, 6.4770e-02,
		6.7246e-02, 6.9258e-02, 7.0800e-02, 7.1872e-02, 7.2481e-02, 7.2639e-02,
		7.2365e-02, 7.1681e-02, 7.0617e-02, 6.9204e-02, 6.7480e-02, 6.5481e-02,
		6.3248e-02, 6.0821e-02, 5.8240e-02, 5.5544e-02, 5.2769e-02, 4.9950e-02,
		4.7118e-02, 4.4302e-02, 4.1525e-02, 3.8810e-02, 3.6173e-02, 3.3630e-02,
		3.1190e-02, 2.8863e-02, 2.6654e-02, 2.4567e-02, 2.2602e-02, 2.0759e-02,
		1.9038e-02, 1.7434e-02, 1.5944e-02, 1.4564e-02, 1.3288e-02, 1.2112e-02,
		1.1030e-02, 1.0036e-02, 9.1242e-03, 8.2899e-03, 7.5272e-03, 6.8310e-03,
		6.1962e-03, 5.6182e-03 };
	
    // ==================================
    // == Thermal photons, event-by-event ideal hydro
    // == Chatterjee, Holopainen, Renk and Eskola, Phys. Rev. C85 (2012) 064910
    // ==================================
    
    const int nPointsChatterjee = 6;
    const double pTChatterjee[nPointsChatterjee] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    
    const double thermalChatterjee00to20[nPointsChatterjee] = {
        0.58458388, 0.034374192, 0.0035716135, 0.00051872584, 9.3080642e-05, 1.9133547e-05
    };
    
    const double thermalChatterjee20to40[nPointsChatterjee] = {
        0.23079249, 0.013742025, 0.0014018224, 0.00019428198, 3.2901498e-05, 6.3852253e-06
    };
    
    const double thermalChatterjee40to60[nPointsChatterjee] = {
        0.079274714, 0.0046140193, 0.00044801761, 5.7851961e-05, 9.0493122e-06, 1.6211998e-06
    };
    
    // ==================================
    // == Thermal and pQCD photons, event-by-event ideal hydro & 
    // == Chatterjee, Holopainen, Helenius, Renk and Eskola, Phys. Rev. C88 (2013) 034901
    // ==================================
    
    const double thermalChatterjee00to40[nPointsChatterjee] = {
        0.36408135, 0.021779899, 0.0022418513, 0.00031302558, 5.2974472e-05, 1.0172486e-05
    };
    
    // NOTE! Be careful, number of points the same but the first pT -bin is different
    //       when NLO pQCD photons included!
    const double pTprimeChatterjee[nPointsChatterjee] = {1.3, 2.0, 3.0, 4.0, 5.0, 6.0};
    
    // here "pQCD" = "prompt + fragmentation" @ NLO, using EPS09s
    const double pQCDChatterjee00to40[nPointsChatterjee] = {
        0.056673277, 0.010903951, 0.0017116971, 0.00042685796, 0.0001460677, 5.8856003e-05
    };
    
    // here "direct" means "thermal + pQCD"
    const double directChatterjee00to40[nPointsChatterjee] = {
        0.21327, 0.032683849, 0.003953548, 0.000739884, 0.000199042, 6.902849e-05
    };
    
    // ==================================
    // == vn's of thermal photons, event-by-event ideal hydro
    // == Chatterjee, Srivastava, Renk, arXiv:1401.7464 [hep-ph]
    // ==================================
    
    
    // PP = Participant Plane, FIC = Fluctuating Initial Conditions
    double v2ChatterjeePPFIC00to40[nPointsChatterjee] = {
        0.019619884, 0.029385347, 0.028213995, 0.022445135, 0.016158229, 0.01124153
    };
    
    // RP = Reaction Plane
    double v2ChatterjeeRPFIC00to40[nPointsChatterjee] = {
        0.010629522, 0.018450325, 0.020289814, 0.017597947, 0.013535189, 0.0096963411
    };

    // SIC = Smooth Initial Conditions
    double v2ChatterjeeSIC00to40[nPointsChatterjee] = {
        0.01942, 0.02251, 0.01217, 0.004691, 0.001656, 0.000643
    };
    
    double v3ChatterjeePPFIC00to40[nPointsChatterjee] = {
        0.00995129067, 0.0204267316, 0.0208320301, 0.0166724697, 0.0120276138, 0.00829880964
    };
    
    // ==================================
    // == pQCD photons (prompt + fragmentation) using EPS09s
    // == Helenius, Eskola, Paukkunen, JHEP 1305 (2013) 030 (<--- calculation itself)
    // == Chatterjee, Holopainen, Helenius, Renk and Eskola, Phys. Rev. C88 (2013) 034901 (<--- with thermal emission)
    // ==================================
    
    const int nPointsHelenius = 21;
    const double pTHelenius[nPointsHelenius] = {
        1.3, 2.0, 3.0, 4.0, 5.0, 6.0, 7.5, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0,
        40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0
    };
    
    // 00 - 20 %
    const double pQCDHelenius00to20[nPointsHelenius] = {
        0.0831538, 0.0160724, 0.00253054, 0.00063178, 0.000216407, 8.72282e-05, 2.92529e-05,
        6.95488e-06, 8.84239e-07, 1.98359e-07, 6.07403e-08, 2.27279e-08, 9.77465e-09, 4.65378e-09,
        2.40004e-09, 1.31728e-09, 4.57895e-10, 1.83691e-10, 8.1956e-11, 3.9689e-11, 2.0515e-11
    };

    const double pQCDHelenius00to20High[nPointsHelenius] = {
        0.108298, 0.0195222, 0.00291263, 0.000711682, 0.000239186, 9.58307e-05, 3.1634e-05,
        7.38036e-06, 9.25639e-07, 2.05946e-07, 6.33679e-08, 2.3843e-08, 1.02865e-08, 4.91406e-09,
        2.53445e-09, 1.39018e-09, 4.83129e-10, 1.9365e-10, 8.6213e-11, 4.16319e-11, 2.1441e-11
    };
    
    const double pQCDHelenius00to20Low[nPointsHelenius] = {
        0.0562462, 0.012959, 0.00213725, 0.000553588, 0.000191324, 7.8984e-05, 2.67933e-05, 6.48788e-06,
        8.48868e-07, 1.90677e-07, 5.83731e-08, 2.17816e-08, 9.33357e-09, 4.43725e-09, 2.28039e-09, 1.24847e-09,
        4.33856e-10, 1.74363e-10, 7.79854e-11, 3.78436e-11, 1.95729e-11
    };
    
    // 20 - 40 %
    const double pQCDHelenius20to40[nPointsHelenius] = {
        0.0302371, 0.00574401, 0.000894192, 0.000222269, 7.58425e-05, 3.05298e-05, 1.02002e-05,
        2.41488e-06, 3.05582e-07, 6.82319e-08, 2.08404e-08, 7.78142e-09, 3.34119e-09, 1.58927e-09,
        8.18668e-10, 4.48851e-10, 1.55836e-10, 6.2488e-11, 2.78751e-11, 1.34971e-11, 6.97621e-12
    };
    
    const double pQCDHelenius20to40High[nPointsHelenius] = {
        0.038682, 0.00687381, 0.00102394, 0.000249278, 8.35966e-05, 3.33872e-05, 1.10008e-05,
        2.56353e-06, 3.19264e-07, 7.09122e-08, 2.17234e-08, 8.14425e-09, 3.50565e-09, 1.67144e-09,
        8.61399e-10, 4.72426e-10, 1.64057e-10, 6.57044e-11, 2.92554e-11, 1.41379e-11, 7.28954e-12
    };
    
    const double pQCDHelenius20to40Low[nPointsHelenius] = {
        0.0214729, 0.00465121, 0.000763821, 0.000195246, 6.75528e-05, 2.76775e-05, 9.37705e-06,
        2.26286e-06, 2.93016e-07, 6.56889e-08, 2.00314e-08, 7.45233e-09, 3.18526e-09, 1.51135e-09,
        7.76267e-10, 4.25048e-10, 1.47697e-10, 5.93596e-11, 2.65564e-11, 1.28957e-11, 6.67832e-12
    };
    
    // 40 - 80 %
    const double pQCDHelenius40to80[nPointsHelenius] = {
        0.00511635, 0.000942834, 0.000144338, 3.55969e-05, 1.20717e-05, 4.84382e-06, 1.60866e-06,
        3.78296e-07, 4.73948e-08, 1.05011e-08, 3.19276e-09, 1.18704e-09, 5.07954e-10, 2.41069e-10,
        1.23932e-10, 6.78407e-11, 2.35076e-11, 9.42157e-12, 4.20418e-12, 2.03664e-12, 1.05325e-12
    };
    
    const double pQCDHelenius40to80High[nPointsHelenius] = {
        0.00641301, 0.00110945, 0.000163554, 3.95869e-05, 1.32298e-05, 5.26207e-06, 1.72696e-06,
        4.0046e-07, 4.941e-08, 1.09121e-08, 3.31971e-09, 1.23905e-09, 5.31645e-10, 2.52829e-10,
        1.30095e-10, 7.12687e-11, 2.47114e-11, 9.88763e-12, 4.40169e-12, 2.12821e-12, 1.09873e-12
    };
    
    const double pQCDHelenius40to80Low[nPointsHelenius] = {
        0.00381674, 0.0007756, 0.000124858, 3.15216e-05, 1.0851e-05, 4.41096e-06, 1.48625e-06,
        3.55707e-07, 4.54627e-08, 1.01353e-08, 3.07314e-09, 1.13908e-09, 4.8546e-10, 2.29757e-10,
        1.17848e-10, 6.44702e-11, 2.23743e-11, 8.98612e-12, 4.019e-12, 1.9525e-12, 1.01247e-12
    };
    
    // ================================================
    // == DATA POINTS FROM DIFFERENT PAPERS FINISHED ==
    // == BELOW ONLY GENERATION OF THE TGRAPH'S      ==
    // ================================================
    
	// ================
	// TGraph's for results in Phys.Rev. C84 (2011) 064903
	// ================
    // 00 - 20 %
	gThermalHolopainen00to20      = new TGraph(nPointsHolopainen);
	gThermalHolopainen00to20Tc200 = new TGraph(nPointsHolopainen);
	gThermalHolopainen00to20Bag   = new TGraph(nPointsHolopainen);
	gv2Holopainen00to20           = new TGraph(nPointsHolopainen);
	gv2Holopainen00to20Tc200      = new TGraph(nPointsHolopainen);
	gv2Holopainen00to20Bag        = new TGraph(nPointsHolopainen);
    // 20 - 40 %
	gThermalHolopainen20to40      = new TGraph(nPointsHolopainen);
	gThermalHolopainen20to40Tc200 = new TGraph(nPointsHolopainen);
	gThermalHolopainen20to40Bag   = new TGraph(nPointsHolopainen);
	gv2Holopainen20to40           = new TGraph(nPointsHolopainen);
	gv2Holopainen20to40Tc200      = new TGraph(nPointsHolopainen);
	gv2Holopainen20to40Bag        = new TGraph(nPointsHolopainen);
    // 00 - 40 %
    gThermalHolopainen00to40      = new TGraph(nPointsHolopainen);
	gThermalHolopainen00to40Tc200 = new TGraph(nPointsHolopainen);
	gThermalHolopainen00to40Bag   = new TGraph(nPointsHolopainen);
	gv2Holopainen00to40           = new TGraph(nPointsHolopainen);
	gv2Holopainen00to40Tc200      = new TGraph(nPointsHolopainen);
	gv2Holopainen00to40Bag        = new TGraph(nPointsHolopainen);
    
	
	for(int counter = 0; counter < nPointsHolopainen; counter++){
        // 00 - 40 %
		gThermalHolopainen00to40->SetPoint(counter, pTHolopainen[counter], thermalHolopainen00to40[counter]);
		gThermalHolopainen00to40Tc200->SetPoint(counter, pTHolopainen[counter], thermalHolopainen00to40Tc200[counter]);
		gThermalHolopainen00to40Bag->SetPoint(counter, pTHolopainen[counter], thermalHolopainen00to40Bag[counter]);
		gv2Holopainen00to40->SetPoint(counter, pTHolopainen[counter], v2Holopainen00to40[counter]);
		gv2Holopainen00to40Tc200->SetPoint(counter, pTHolopainen[counter], v2Holopainen00to40Tc200[counter]);
		gv2Holopainen00to40Bag->SetPoint(counter, pTHolopainen[counter], v2Holopainen00to40Bag[counter]);
        // 00 - 20 %
		gThermalHolopainen00to20->SetPoint(counter, pTHolopainen[counter], thermalHolopainen00to20[counter]);
		gThermalHolopainen00to20Tc200->SetPoint(counter, pTHolopainen[counter], thermalHolopainen00to20Tc200[counter]);
		gThermalHolopainen00to20Bag->SetPoint(counter, pTHolopainen[counter], thermalHolopainen00to20Bag[counter]);
		gv2Holopainen00to20->SetPoint(counter, pTHolopainen[counter], v2Holopainen00to20[counter]);
		gv2Holopainen00to20Tc200->SetPoint(counter, pTHolopainen[counter], v2Holopainen00to20Tc200[counter]);
		gv2Holopainen00to20Bag->SetPoint(counter, pTHolopainen[counter], v2Holopainen00to20Bag[counter]);
        // 20 - 40 %
		gThermalHolopainen20to40->SetPoint(counter, pTHolopainen[counter], thermalHolopainen20to40[counter]);
		gThermalHolopainen20to40Tc200->SetPoint(counter, pTHolopainen[counter], thermalHolopainen20to40Tc200[counter]);
		gThermalHolopainen20to40Bag->SetPoint(counter, pTHolopainen[counter], thermalHolopainen20to40Bag[counter]);
		gv2Holopainen20to40->SetPoint(counter, pTHolopainen[counter], v2Holopainen20to40[counter]);
		gv2Holopainen20to40Tc200->SetPoint(counter, pTHolopainen[counter], v2Holopainen20to40Tc200[counter]);
		gv2Holopainen20to40Bag->SetPoint(counter, pTHolopainen[counter], v2Holopainen20to40Bag[counter]);
	}
	
	// Make default plotting settings, print out for check's
    
	// pT-spectra
    
    // Default configurations
	gThermalHolopainen00to40->SetLineStyle(1);
	gThermalHolopainen00to40->SetLineColor(1);
	gThermalHolopainen00to40->SetLineWidth(3);
	gThermalHolopainen00to40->SetTitle("Default 00 - 40 % : EoS L, Tswitch = 170 MeV and TRG rates");
	gThermalHolopainen00to40->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gThermalHolopainen00to40->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
	gThermalHolopainen00to40->Print();
    
	gThermalHolopainen00to20->SetLineStyle(1);
	gThermalHolopainen00to20->SetLineColor(1);
	gThermalHolopainen00to20->SetLineWidth(3);
	gThermalHolopainen00to20->SetTitle("Default 00 - 20 % : EoS L, Tswitch = 170 MeV and TRG rates");
	gThermalHolopainen00to20->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gThermalHolopainen00to20->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
	gThermalHolopainen00to20->Print();
    
	gThermalHolopainen20to40->SetLineStyle(1);
	gThermalHolopainen20to40->SetLineColor(1);
	gThermalHolopainen20to40->SetLineWidth(3);
	gThermalHolopainen20to40->SetTitle("Default 20 - 40 % : EoS L, Tswitch = 170 MeV and TRG rates");
	gThermalHolopainen20to40->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gThermalHolopainen20to40->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
	gThermalHolopainen20to40->Print();
    
    // Larger switching temperature
	gThermalHolopainen00to40Tc200->SetLineStyle(1);
	gThermalHolopainen00to40Tc200->SetLineColor(2);
	gThermalHolopainen00to40Tc200->SetLineWidth(3);
	gThermalHolopainen00to40Tc200->SetTitle("00 - 40 % : EoS L, Tswitch = 200 MeV and TRG rates");
	gThermalHolopainen00to40Tc200->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gThermalHolopainen00to40Tc200->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
	gThermalHolopainen00to40Tc200->Print();
    
	gThermalHolopainen00to20Tc200->SetLineStyle(1);
	gThermalHolopainen00to20Tc200->SetLineColor(2);
	gThermalHolopainen00to20Tc200->SetLineWidth(3);
	gThermalHolopainen00to20Tc200->SetTitle("00 - 20 % : EoS L, Tswitch = 200 MeV and TRG rates");
	gThermalHolopainen00to20Tc200->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gThermalHolopainen00to20Tc200->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
	gThermalHolopainen00to20Tc200->Print();
    
	gThermalHolopainen20to40Tc200->SetLineStyle(1);
	gThermalHolopainen20to40Tc200->SetLineColor(2);
	gThermalHolopainen20to40Tc200->SetLineWidth(3);
	gThermalHolopainen20to40Tc200->SetTitle("20 - 40 % : EoS L, Tswitch = 200 MeV and TRG rates");
	gThermalHolopainen20to40Tc200->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gThermalHolopainen20to40Tc200->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
	gThermalHolopainen20to40Tc200->Print();
    
    // Bag EoS and old rates
	gThermalHolopainen00to40Bag->SetLineStyle(1);
	gThermalHolopainen00to40Bag->SetLineColor(4);
	gThermalHolopainen00to40Bag->SetLineWidth(3);
	gThermalHolopainen00to40Bag->SetTitle("00 - 40 % : EoS Q, Tc = 165 MeV and R92 rates");
	gThermalHolopainen00to40Bag->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gThermalHolopainen00to40Bag->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
	gThermalHolopainen00to40Bag->Print();
    
	gThermalHolopainen00to20Bag->SetLineStyle(1);
	gThermalHolopainen00to20Bag->SetLineColor(4);
	gThermalHolopainen00to20Bag->SetLineWidth(3);
	gThermalHolopainen00to20Bag->SetTitle("00 - 20 % : EoS Q, Tc = 165 MeV and R92 rates");
	gThermalHolopainen00to20Bag->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gThermalHolopainen00to20Bag->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
	gThermalHolopainen00to20Bag->Print();
    
	gThermalHolopainen20to40Bag->SetLineStyle(1);
	gThermalHolopainen20to40Bag->SetLineColor(4);
	gThermalHolopainen20to40Bag->SetLineWidth(3);
	gThermalHolopainen20to40Bag->SetTitle("20 - 40 % : EoS Q, Tc = 165 MeV and R92 rates");
	gThermalHolopainen20to40Bag->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gThermalHolopainen20to40Bag->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
	gThermalHolopainen20to40Bag->Print();
    
	// .. and the same for v2
    
    // Default configuration
	gv2Holopainen00to40->SetLineStyle(1);
	gv2Holopainen00to40->SetLineColor(1);
	gv2Holopainen00to40->SetLineWidth(3);
	gv2Holopainen00to40->SetTitle("Default 00 - 40 % : EoS L, Tswitch = 170 MeV and TRG rates");
	gv2Holopainen00to40->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gv2Holopainen00to40->GetYaxis()->SetTitle("v_{2}");
	gv2Holopainen00to40->Print();
    
	gv2Holopainen00to20->SetLineStyle(1);
	gv2Holopainen00to20->SetLineColor(1);
	gv2Holopainen00to20->SetLineWidth(3);
	gv2Holopainen00to20->SetTitle("Default 00 - 20 % : EoS L, Tswitch = 170 MeV and TRG rates");
	gv2Holopainen00to20->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gv2Holopainen00to20->GetYaxis()->SetTitle("v_{2}");
	gv2Holopainen00to20->Print();
    
	gv2Holopainen20to40->SetLineStyle(1);
	gv2Holopainen20to40->SetLineColor(1);
	gv2Holopainen20to40->SetLineWidth(3);
	gv2Holopainen20to40->SetTitle("Default 20 - 40 % : EoS L, Tswitch = 170 MeV and TRG rates");
	gv2Holopainen20to40->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gv2Holopainen20to40->GetYaxis()->SetTitle("v_{2}");
	gv2Holopainen20to40->Print();
    
    // Larger switching temperature
	gv2Holopainen00to40Tc200->SetLineStyle(1);
	gv2Holopainen00to40Tc200->SetLineColor(2);
	gv2Holopainen00to40Tc200->SetLineWidth(3);
	gv2Holopainen00to40Tc200->SetTitle("00 - 40 % : EoS L, Tswitch = 200 MeV and TRG rates");
	gv2Holopainen00to40Tc200->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gv2Holopainen00to40Tc200->GetYaxis()->SetTitle("v_{2}");
	gv2Holopainen00to40Tc200->Print();
    
	gv2Holopainen00to20Tc200->SetLineStyle(1);
	gv2Holopainen00to20Tc200->SetLineColor(2);
	gv2Holopainen00to20Tc200->SetLineWidth(3);
	gv2Holopainen00to20Tc200->SetTitle("00 - 20 % : EoS L, Tswitch = 200 MeV and TRG rates");
	gv2Holopainen00to20Tc200->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gv2Holopainen00to20Tc200->GetYaxis()->SetTitle("v_{2}");
	gv2Holopainen00to20Tc200->Print();
    
	gv2Holopainen20to40Tc200->SetLineStyle(1);
	gv2Holopainen20to40Tc200->SetLineColor(2);
	gv2Holopainen20to40Tc200->SetLineWidth(3);
	gv2Holopainen20to40Tc200->SetTitle("20 - 40 % : EoS L, Tswitch = 200 MeV and TRG rates");
	gv2Holopainen20to40Tc200->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gv2Holopainen20to40Tc200->GetYaxis()->SetTitle("v_{2}");
	gv2Holopainen20to40Tc200->Print();
    
    // Bag EoS and old rates
	gv2Holopainen00to40Bag->SetLineStyle(1);
	gv2Holopainen00to40Bag->SetLineColor(4);
	gv2Holopainen00to40Bag->SetLineWidth(3);
	gv2Holopainen00to40Bag->SetTitle("00 - 40 % : EoS Q, Tc = 165 MeV and R92 rates");
	gv2Holopainen00to40Bag->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gv2Holopainen00to40Bag->GetYaxis()->SetTitle("v_{2}");
    gv2Holopainen00to40Bag->Print();
	
	gv2Holopainen00to20Bag->SetLineStyle(1);
	gv2Holopainen00to20Bag->SetLineColor(4);
	gv2Holopainen00to20Bag->SetLineWidth(3);
	gv2Holopainen00to20Bag->SetTitle("00 - 20 % : EoS Q, Tc = 165 MeV and R92 rates");
	gv2Holopainen00to20Bag->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gv2Holopainen00to20Bag->GetYaxis()->SetTitle("v_{2}");
    gv2Holopainen00to20Bag->Print();
	
	gv2Holopainen20to40Bag->SetLineStyle(1);
	gv2Holopainen20to40Bag->SetLineColor(4);
	gv2Holopainen20to40Bag->SetLineWidth(3);
	gv2Holopainen20to40Bag->SetTitle("20 - 40 % : EoS Q, Tc = 165 MeV and R92 rates");
	gv2Holopainen20to40Bag->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gv2Holopainen20to40Bag->GetYaxis()->SetTitle("v_{2}");
    gv2Holopainen20to40Bag->Print();
	
	// ================
	// TGraph's for results in papers:
    // Phys. Rev. C85 (2012) 064910
    // Phys. Rev. C88 (2013) 034901
    // arXiv:1401.7464 [hep-ph]
	// ================
    
    // spectra
    gThermalChatterjee00to20 = new TGraph(nPointsChatterjee);
    gThermalChatterjee20to40 = new TGraph(nPointsChatterjee);
    gThermalChatterjee40to60 = new TGraph(nPointsChatterjee);
    gThermalChatterjee00to40 = new TGraph(nPointsChatterjee);
    gpQCDChatterjee00to40    = new TGraph(nPointsChatterjee);
    gDirectChatterjee00to40  = new TGraph(nPointsChatterjee);
    
    // v2 and v3
    gv2ChatterjeePPFIC00to40 = new TGraph(nPointsChatterjee);
    gv2ChatterjeeRPFIC00to40 = new TGraph(nPointsChatterjee);
    gv2ChatterjeeSIC00to40   = new TGraph(nPointsChatterjee);
    gv3ChatterjeePPFIC00to40 = new TGraph(nPointsChatterjee);
    
    for(int counter = 0; counter < nPointsChatterjee; counter++){
        double pT = pTChatterjee[counter];
        
        // spectra
        gThermalChatterjee00to20->SetPoint(counter, pT, thermalChatterjee00to20[counter]);
        gThermalChatterjee20to40->SetPoint(counter, pT, thermalChatterjee20to40[counter]);
        gThermalChatterjee40to60->SetPoint(counter, pT, thermalChatterjee40to60[counter]);
        gThermalChatterjee00to40->SetPoint(counter, pT, thermalChatterjee00to40[counter]);
        // NOTE: different pT's ! In purpose.
        gpQCDChatterjee00to40->SetPoint(counter, pTprimeChatterjee[counter], pQCDChatterjee00to40[counter]);
        gDirectChatterjee00to40->SetPoint(counter, pTprimeChatterjee[counter], directChatterjee00to40[counter]);
        
        // v2 and v3
        gv2ChatterjeePPFIC00to40->SetPoint(counter, pT, v2ChatterjeePPFIC00to40[counter]);
        gv2ChatterjeeRPFIC00to40->SetPoint(counter, pT, v2ChatterjeeRPFIC00to40[counter]);
        gv2ChatterjeeSIC00to40->SetPoint(counter, pT, v2ChatterjeeSIC00to40[counter]);
        gv3ChatterjeePPFIC00to40->SetPoint(counter, pT, v3ChatterjeePPFIC00to40[counter]);
    }
    
    gThermalChatterjee00to20->SetLineStyle(1);
    gThermalChatterjee00to20->SetLineColor(1);
    gThermalChatterjee00to20->SetMarkerStyle(20);
    gThermalChatterjee00to20->SetMarkerSize(1.5);
    gThermalChatterjee00to20->SetMarkerColor(1);
    gThermalChatterjee00to20->SetLineWidth(1);
    gThermalChatterjee00to20->SetTitle("Thermal photons - centrality 0 - 20 %");
	gThermalChatterjee00to20->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gThermalChatterjee00to20->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
    gThermalChatterjee00to20->Print();
    
    gThermalChatterjee20to40->SetLineStyle(1);
    gThermalChatterjee20to40->SetLineColor(2);
    gThermalChatterjee20to40->SetMarkerStyle(21);
    gThermalChatterjee20to40->SetMarkerSize(1.5);
    gThermalChatterjee20to40->SetMarkerColor(2);
    gThermalChatterjee20to40->SetLineWidth(1);
    gThermalChatterjee20to40->SetTitle("Thermal photons - centrality 20 - 40 %");
	gThermalChatterjee20to40->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gThermalChatterjee20to40->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
    gThermalChatterjee20to40->Print();
    
    gThermalChatterjee40to60->SetLineStyle(1);
    gThermalChatterjee40to60->SetLineColor(6);
    gThermalChatterjee40to60->SetMarkerStyle(34);
    gThermalChatterjee40to60->SetMarkerSize(1.5);
    gThermalChatterjee40to60->SetMarkerColor(6);
    gThermalChatterjee40to60->SetLineWidth(1);
    gThermalChatterjee40to60->SetTitle("Thermal photons - centrality 40 - 60 %");
	gThermalChatterjee40to60->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gThermalChatterjee40to60->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
    gThermalChatterjee40to60->Print();
    
    gThermalChatterjee00to40->SetLineStyle(1);
    gThermalChatterjee00to40->SetLineColor(4);
    gThermalChatterjee00to40->SetMarkerStyle(33);
    gThermalChatterjee00to40->SetMarkerSize(1.5);
    gThermalChatterjee00to40->SetMarkerColor(4);
    gThermalChatterjee00to40->SetLineWidth(1);
    gThermalChatterjee00to40->SetTitle("Thermal photons - centrality 00 - 40 %");
	gThermalChatterjee00to40->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gThermalChatterjee00to40->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
    gThermalChatterjee00to40->Print();
    
    gpQCDChatterjee00to40->SetLineStyle(2);
    gpQCDChatterjee00to40->SetLineColor(4);
    gpQCDChatterjee00to40->SetMarkerStyle(28);
    gpQCDChatterjee00to40->SetMarkerSize(1.5);
    gpQCDChatterjee00to40->SetMarkerColor(4);
    gpQCDChatterjee00to40->SetLineWidth(1);
    gpQCDChatterjee00to40->SetTitle("pQCD photons - centrality 00 - 40 %");
	gpQCDChatterjee00to40->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gpQCDChatterjee00to40->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
    gpQCDChatterjee00to40->Print();
    
    gDirectChatterjee00to40->SetLineStyle(1);
    gDirectChatterjee00to40->SetLineColor(1);
    gDirectChatterjee00to40->SetMarkerStyle(34);
    gDirectChatterjee00to40->SetMarkerSize(1.5);
    gDirectChatterjee00to40->SetMarkerColor(1);
    gDirectChatterjee00to40->SetLineWidth(1);
    gDirectChatterjee00to40->SetTitle("Direct photons - centrality 00 - 40 %");
	gDirectChatterjee00to40->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gDirectChatterjee00to40->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
    gDirectChatterjee00to40->Print();
    
    // v2 and v3
    gv2ChatterjeePPFIC00to40->SetLineStyle(1);
    gv2ChatterjeePPFIC00to40->SetLineColor(1);
    gv2ChatterjeePPFIC00to40->SetMarkerStyle(20);
    gv2ChatterjeePPFIC00to40->SetMarkerSize(1.5);
    gv2ChatterjeePPFIC00to40->SetMarkerColor(1);
    gv2ChatterjeePPFIC00to40->SetLineWidth(1);
    gv2ChatterjeePPFIC00to40->SetTitle("Thermal photon v2 - centrality 0 - 40 %, wrt PP plane");
	gv2ChatterjeePPFIC00to40->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gv2ChatterjeePPFIC00to40->GetYaxis()->SetTitle("thermal photon v_{2}(PP)");
    gv2ChatterjeePPFIC00to40->Print();
    
    gv2ChatterjeeRPFIC00to40->SetLineStyle(1);
    gv2ChatterjeeRPFIC00to40->SetLineColor(2);
    gv2ChatterjeeRPFIC00to40->SetMarkerStyle(21);
    gv2ChatterjeeRPFIC00to40->SetMarkerSize(1.5);
    gv2ChatterjeeRPFIC00to40->SetMarkerColor(2);
    gv2ChatterjeeRPFIC00to40->SetLineWidth(1);
    gv2ChatterjeeRPFIC00to40->SetTitle("Thermal photon v2 - centrality 0 - 40 %, wrt RP plane");
	gv2ChatterjeeRPFIC00to40->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gv2ChatterjeeRPFIC00to40->GetYaxis()->SetTitle("thermal photon v_{2}(RP)");
    gv2ChatterjeeRPFIC00to40->Print();
    
    gv2ChatterjeeSIC00to40->SetLineStyle(1);
    gv2ChatterjeeSIC00to40->SetLineColor(4);
    gv2ChatterjeeSIC00to40->SetMarkerStyle(33);
    gv2ChatterjeeSIC00to40->SetMarkerSize(1.5);
    gv2ChatterjeeSIC00to40->SetMarkerColor(4);
    gv2ChatterjeeSIC00to40->SetLineWidth(1);
    gv2ChatterjeeSIC00to40->SetTitle("Thermal photon v2 - centrality 0 - 40 %, SIC");
	gv2ChatterjeeSIC00to40->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gv2ChatterjeeSIC00to40->GetYaxis()->SetTitle("thermal photon v_{2}");
    gv2ChatterjeeSIC00to40->Print();
    
    
    gv3ChatterjeePPFIC00to40->SetLineStyle(1);
    gv3ChatterjeePPFIC00to40->SetLineColor(1);
    gv3ChatterjeePPFIC00to40->SetMarkerStyle(20);
    gv3ChatterjeePPFIC00to40->SetMarkerSize(1.5);
    gv3ChatterjeePPFIC00to40->SetMarkerColor(1);
    gv3ChatterjeePPFIC00to40->SetLineWidth(1);
    gv3ChatterjeePPFIC00to40->SetTitle("Thermal photon v3 - centrality 0 - 40 %, wrt PP plane");
	gv3ChatterjeePPFIC00to40->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gv3ChatterjeePPFIC00to40->GetYaxis()->SetTitle("thermal photon v_{3}(PP)");
    gv3ChatterjeePPFIC00to40->Print();

    // ================
    // TGraph's for results in JHEP 1305 (2013) 030
    // performed for centrality bins asked
    // ================
    
    // EPS09s central set
    gpQCDHelenius00to20     = new TGraph(nPointsHelenius);
    gpQCDHelenius20to40     = new TGraph(nPointsHelenius);
    gpQCDHelenius40to80     = new TGraph(nPointsHelenius);
    // EPS09s error sets - high limit
    gpQCDHeleniusHigh00to20 = new TGraph(nPointsHelenius);
    gpQCDHeleniusHigh20to40 = new TGraph(nPointsHelenius);
    gpQCDHeleniusHigh40to80 = new TGraph(nPointsHelenius);
    // EPS09s error sets - low limit
    gpQCDHeleniusLow00to20  = new TGraph(nPointsHelenius);
    gpQCDHeleniusLow20to40  = new TGraph(nPointsHelenius);
    gpQCDHeleniusLow40to80  = new TGraph(nPointsHelenius);
    
    for (int counter = 0; counter < nPointsHelenius; counter++){
        double pT = pTHelenius[counter];
        // central
        gpQCDHelenius00to20->SetPoint(counter, pT, pQCDHelenius00to20[counter]);
        gpQCDHelenius20to40->SetPoint(counter, pT, pQCDHelenius20to40[counter]);
        gpQCDHelenius40to80->SetPoint(counter, pT, pQCDHelenius40to80[counter]);
        // high
        gpQCDHeleniusHigh00to20->SetPoint(counter, pT, pQCDHelenius00to20High[counter]);
        gpQCDHeleniusHigh20to40->SetPoint(counter, pT, pQCDHelenius20to40High[counter]);
        gpQCDHeleniusHigh40to80->SetPoint(counter, pT, pQCDHelenius40to80High[counter]);
        // low
        gpQCDHeleniusLow00to20->SetPoint(counter, pT, pQCDHelenius00to20Low[counter]);
        gpQCDHeleniusLow20to40->SetPoint(counter, pT, pQCDHelenius20to40Low[counter]);
        gpQCDHeleniusLow40to80->SetPoint(counter, pT, pQCDHelenius40to80Low[counter]);
    }
    
    gpQCDHelenius00to20->SetLineStyle(1);
    gpQCDHelenius00to20->SetLineColor(1);
    gpQCDHelenius00to20->SetLineWidth(3);
    gpQCDHelenius00to20->SetTitle("centrality 0 - 20 % - central EPS09s set");
	gpQCDHelenius00to20->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gpQCDHelenius00to20->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
    gpQCDHelenius00to20->Print();
    
    gpQCDHeleniusHigh00to20->SetLineStyle(2);
    gpQCDHeleniusHigh00to20->SetLineColor(1);
    gpQCDHeleniusHigh00to20->SetLineWidth(1);
    gpQCDHeleniusHigh00to20->SetTitle("centrality 0 - 20 % - high EPS09s set");
	gpQCDHeleniusHigh00to20->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gpQCDHeleniusHigh00to20->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
    gpQCDHeleniusHigh00to20->Print();
    
    gpQCDHeleniusLow00to20->SetLineStyle(2);
    gpQCDHeleniusLow00to20->SetLineColor(1);
    gpQCDHeleniusLow00to20->SetLineWidth(1);
    gpQCDHeleniusLow00to20->SetTitle("centrality 0 - 20 % - low EPSs set");
	gpQCDHeleniusLow00to20->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gpQCDHeleniusLow00to20->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
    gpQCDHeleniusLow00to20->Print();
    
    gpQCDHelenius20to40->SetLineStyle(1);
    gpQCDHelenius20to40->SetLineColor(2);
    gpQCDHelenius20to40->SetLineWidth(3);
    gpQCDHelenius20to40->SetTitle("centrality 20 - 40 % - central EPS09s set");
	gpQCDHelenius20to40->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gpQCDHelenius20to40->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
    gpQCDHelenius20to40->Print();
    
    gpQCDHeleniusHigh20to40->SetLineStyle(2);
    gpQCDHeleniusHigh20to40->SetLineColor(2);
    gpQCDHeleniusHigh20to40->SetLineWidth(1);
    gpQCDHeleniusHigh20to40->SetTitle("centrality 20 - 40 % - high EPS09s set");
	gpQCDHeleniusHigh20to40->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gpQCDHeleniusHigh20to40->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
    gpQCDHeleniusHigh20to40->Print();
    
    gpQCDHeleniusLow20to40->SetLineStyle(2);
    gpQCDHeleniusLow20to40->SetLineColor(2);
    gpQCDHeleniusLow20to40->SetLineWidth(1);
    gpQCDHeleniusLow20to40->SetTitle("centrality 20 - 40 % - low EPSs set");
	gpQCDHeleniusLow20to40->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gpQCDHeleniusLow20to40->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
    gpQCDHeleniusLow20to40->Print();
    
    gpQCDHelenius40to80->SetLineStyle(1);
    gpQCDHelenius40to80->SetLineColor(4);
    gpQCDHelenius40to80->SetLineWidth(3);
    gpQCDHelenius40to80->SetTitle("centrality 40 - 80 % - central EPS09s set");
	gpQCDHelenius40to80->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gpQCDHelenius40to80->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
    gpQCDHelenius40to80->Print();
    
    gpQCDHeleniusHigh40to80->SetLineStyle(2);
    gpQCDHeleniusHigh40to80->SetLineColor(4);
    gpQCDHeleniusHigh40to80->SetLineWidth(1);
    gpQCDHeleniusHigh40to80->SetTitle("centrality 40 - 80 % - high EPS09s set");
	gpQCDHeleniusHigh40to80->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gpQCDHeleniusHigh40to80->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
    gpQCDHeleniusHigh40to80->Print();
    
    gpQCDHeleniusLow40to80->SetLineStyle(2);
    gpQCDHeleniusLow40to80->SetLineColor(4);
    gpQCDHeleniusLow40to80->SetLineWidth(1);
    gpQCDHeleniusLow40to80->SetTitle("centrality 40 - 80 % - low EPSs set");
	gpQCDHeleniusLow40to80->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	gpQCDHeleniusLow40to80->GetYaxis()->SetTitle("E dN/d^{3}p at y = 0  [(GeV/c)^{-2}]");
    gpQCDHeleniusLow40to80->Print();

}