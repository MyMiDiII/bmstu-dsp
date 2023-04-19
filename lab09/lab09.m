function main()
I=double(imread('bimage2.bmp')) / 255;

figure;
imshow(I);

PSF=fspecial('motion', 50, 25);
J1=deconvblind(I, PSF);

figure;
imshow(J1);
end
