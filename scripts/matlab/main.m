ppm_list = {'000', '001' '002', '003', '004', '005', '006', '007', '008', '009', '010', '011', '012'};
Y_PNSR_SUM = 0;
YUV_PNSR_SUM = 0;
Y_SSIM_SUM = 0;

for i = 1:numel(ppm_list)
    for j = 1:numel(ppm_list)
        arq1 = ['/home/douglascorrea/light-field/dataset/Lenslet/Bikes/', ppm_list{i}, '_', ppm_list{j},'.ppm'];
        arq2 = ['/home/douglascorrea/light-field/resultados_prediction/Bikes/decoding/', ppm_list{i}, '_', ppm_list{j},'.ppm'];
        disp(arq1);
        disp(arq2);
        im1=imread(arq1);
        im2=imread(arq2);

        ref=bitshift(im1,-6); rec=bitshift(im2,-6);

        [Y_PSNR YUV_PSNR Y_SSIM] = QM(ref,rec,10,10);

        Y_PNSR_SUM = Y_PNSR_SUM + Y_PSNR;
        YUV_PNSR_SUM = YUV_PNSR_SUM + YUV_PSNR;
        Y_SSIM_SUM = Y_SSIM_SUM + Y_SSIM;
    end;
end;

AVERAGE_Y_PNSR = Y_PNSR_SUM/169;
AVERAGE_YUV_PNSR =  YUV_PNSR_SUM/169;
AVERAGE_Y_SSIM =  Y_SSIM_SUM/169;