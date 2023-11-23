//
// Created by JohnLin on 2023/11/17.
//

#ifndef AWGN_H
#define AWGN_H

double get_awgn_noise(double sigma);
double get_sigma_with_snr(double snr_db, double code_rate);

#endif // AWGN_H
