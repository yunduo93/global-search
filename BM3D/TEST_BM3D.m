clc;
clear;
Image=imread('y_est.bmp');
[PSNR, y_est] = CBM3D(1, im2double(Image), 5 , 'np', 1);