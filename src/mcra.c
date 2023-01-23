#include <math.h>
#include <stdio.h>
#include <string.h>
#include "define.h"
#include "DSP_denoise.h"
#include "Platform.h"
#include "denoise_function.h"
#include "table.h"
#include "fft.h"


void mcra_init(){

for (int i = 1; i < (WIN_ALL + 1); ++i){   //
	hann_win[i - 1] = (0.5 - 0.5 * cos(2.0 * PI*(i) / (WIN_ALL + 1)));// 0-----1
	//hann_win[i] = sin(.5*M_PI*sin(.5*M_PI*(i+.5)/FRAME_SIZE) * sin(.5*M_PI*(i+.5)/FRAME_SIZE));
}


for (int i = 0; i <= C_F; i++) {
	pr_snr[i] = PR_SNR;  //可以取更大值？？？ 参考matlab

}
}
static int expintpow_solution(float v_subscript) {
	int vec = 0;
	int g = 0;
	v_subscript *= 100;
	vec = ((int)v_subscript);  // / 0.0001;
	vec = vec < 1 ? 1 : vec;
	vec = vec > 1500 ? 1500 : vec;
	 
	//g = (m_int_value[vec - 1]) ;
	g= (int_value[vec - 1]);
	return g;
}

static int subexp_solution(float v_subscript) {

	int vec = 0;
	int g = 0;
	v_subscript *= 100;
	vec = ((int)v_subscript);  // / 0.0001;
	vec = vec < 1 ? 1 : vec;
	vec = vec > 1500 ? 1500 : vec;
	g = (expsub_value[vec - 1]);
	return g;
}

void NoiseEstimation(int Block){
    int p;
    int L;
    if (Block >= 2) {
        L = 1;

    }
	else {
		L = 1;
		if (Block == 0) {
			for (int i = 0; i <= (C_F); i++) {   //对4000hz以下 估计前导噪声
				amp[i] = amp_pr[i];
				amp_min[i] = amp[i];
				amp_tmp[i] = amp[i];
				noise_est[i] = amp[i];
			}
			memset(init_p, 0, sizeof(int)*((C_F)+1)); //
		}
	}

    //更新估计的噪声
	if ((Block + 1) % L == 0) {
		for (int k = 0; k <= ((C_F)); k++) {
			amp_min[k] = MIN(amp_tmp[k], amp[k]);
			amp_tmp[k] = amp[k];
		}
	}
	
    
    for (int k = 0; k <= ((C_F)); k++)
	{
		//对估计的噪声进行平滑处理
        amp[k] = (0.8*amp[k] + 0.2*amp_pr[k] );//
		amp_min[k] = MIN(amp_min[k], amp[k]);
		amp_tmp[k] = MIN(amp_tmp[k], amp[k]);

        //二次平滑， P是条件语音存在概率
		if (amp[k] > (3*amp_min[k])) //20~~~~~~21
			p = 1;
		else
			p = 0;
        init_p[k] = (0.4 * init_p[k]  + 0.6 * p); //21
		noise_est[k] = ((0.9*noise_est[k] + 0.1*amp_pr[k])*(1 - init_p[k]) + noise_est[k] * init_p[k]); //18


        // MCRA部分 并没有用到信噪比相关的参数 也没有信噪比 * 窗函数操作
	}    
}



void SpeechAbsenceEstm(){

    float sum = 0, a = 0;
    float p_frame, mu, snr_peak;
    float old_snr_frame = 0, snr_frame = 0;

    // 先验信噪比时域平滑  在后面用
	for (int k = 0; k <= (C_F); k++) {   //
		snr[k] = (0.75 * old_snr[k] + 0.25 * current_snr[k]);
	}

    for (int k = 0; k <= (C_F); k++) {

        //前WIN_HALF个频点默认是全局信噪比了
		if (k <= WIN_HALF - 1) {  // 对一帧数据掐头去尾
			snr_global[k] = snr[k];

            //前WIN_HALF个频点 局部单频点信噪比依然用了平滑了 ----3值平滑
			if (k == 0)
				snr_local[k] = snr[k];
			else
				snr_local[k] = (snr[k - 1] + 2 * snr[k] + snr[k + 1]) / 4.f;
		}
		else if (k >= ((C_F) - (WIN_HALF+1))){
				snr_local[k] = snr[k];
				snr_local[k] = (snr[k - 1] + 2 * snr[k] + snr[k + 1]) / 4.f;
		}
		else {
			float a = 0;
            //局部信噪比 窗长
			for (int j = 0; j < WIN_ALL; j++) //点 ***窗函数 反过来相乘了
				a = a + ((hann_win[j] * snr[k + WIN_HALF - j]));
			snr_global[k] = ((a / (WIN_HALF + 1)));

            //依然是3频点平滑
			snr_local[k] = (snr[k - 1] + 2 * snr[k] + snr[k + 1]) / 4.f;
		}

        //语音存在可能性 局部判定
		if (snr_local[k] <= SNR_MIN)   //% (14)
			p_local[k] = 0;
		else if (snr_local[k] >= SNR_MAX)
			p_local[k] = 1;
		else
			p_local[k] = 1 * (log(snr_local[k]) - log(SNR_MIN)) / (log(SNR_MAX)-log(SNR_MIN));

        // 语音存在可能性 全局判定
		if (snr_global[k] <= SNR_MIN)   //% (14)
			p_global[k] = 0;
		else if (snr_global[k] >= SNR_MAX)
			p_global[k] = 1;
		else
			p_global[k] = 1 * (log(snr_global[k]) - log(SNR_MIN)) / (log(SNR_MAX)-log(SNR_MIN));

        //信噪比加起来 后面要用
		sum += snr[k];
	}
    snr_frame = sum / ((C_F)+1);
	sum = 0 ;
	snr_peak = MIN(MAX(snr_frame, SNR_P_MIN ), SNR_P_MAX );//10000

	if (snr_frame <= (snr_peak * SNR_MIN))   // (15)
		mu = 0;
	else if (snr_frame >= (snr_peak * SNR_MAX))
		mu = 1;
	else
		mu = 1 * log(snr_frame  * 1 / snr_peak / SNR_MIN) / log(SNR_MAX / SNR_MIN);

	if (snr_frame > SNR_MIN)
	{
		if (snr_frame > old_snr_frame)
			p_frame = 1;
		else
		{
			p_frame = mu;
		}
    }
    else
	{
		p_frame = 0;
	}
	for (int k = 0; k <= (C_F); k++)
	{
		q_est[k] = 1 - (p_local[k] * p_global[k] * p_frame );//10000  (16)
		q_est[k] = MIN(q_est[k], 0.92);
		old_snr[k] = snr[k];
		//m_old_cosen[FRAME_SIZE - k] = m_old_cosen[k];
	}
    
}

float G_calculate_process(DenoiseState *st, complex *X, complex *DSP_X, int Block) {  // -49

    //frame_analysis(st, X, in);

	int post_temp;
	for (int i = 0; i <= (C_F); i++) {
		amp_pr[i] = sqrt(pow(X[i].r, 2) + pow(X[i].i, 2));  //
		//_EN_cos[i] = ((__int64)winData[i].real << 12) / (1 > m_abs_Y[i] ? 1 : m_abs_Y[i]);
		//m_EN_sin[i] = ((__int64)winData[i].imag << 12) / (1 > m_abs_Y[i] ? 1 : m_abs_Y[i]);
	}
	mcra_init();
    //这里求出的是全频带的噪声 
	NoiseEstimation(Block); 
	for (int i = 0; i <= (C_F); i++) {
        //先求出后验信噪比
		post_snr[i] = MIN((pow(((amp_pr[i]) / noise_est[i]), 2)), 10000); //40db
		
        post_temp = MAX(post_snr[i] - 1, 0);
		
        //DD准则 求先验信噪比 全频带的
        current_snr[i] = MIN(MAX((0.8 * pr_snr[i] + 0.2 * post_temp ), P_MIN), 10000); ///40db
		//m_E_pr_SNR[FRAME_SIZE - i] = m_E_pr_SNR[i];                                                // 0.0001 * (2^24=16777216) =1678  167772 
		
        //有先验信噪比 和后验信噪比 根据查表法求出 GH1 这个是语音存在时候的概率
        v_int[i] = MAX(MIN((current_snr[i] * post_snr[i] ) / (1+ current_snr[i]), 10000), 0.01); // (7~~~8)
		integra[i] = expintpow_solution(v_int[i]);  //
		m_int[i] = integra[i] /16384.f ;
		//m_int[i] = 0.5 * expint(m_v[i])
	

		gh1[i] = MIN(current_snr[i] * m_int[i] / (1+ current_snr[i]), 8);  //(7)
	}

    //语音不存在概率
	SpeechAbsenceEstm();
	for (int i = 0; i <= (FRAME_SIZE); i++) {
        //指数查表
		if(0<i && i<C_F){
		integra[i] = subexp_solution(v_int[i]);  //  (11)
		m_int[i] = integra[i] /16384.f ;
		//m_int[i] = exp(-v[i]);
		arr_temp[i] = 1 + ((1 + current_snr[i]) * integra[i] * q_est[i])  / (1- q_est[i]);

		
        p_est[i] = MAX(1 / arr_temp[i], 0.0001);
		g[i] = pow(gh1[i], p_est[i] ) * pow(NOISE_FACTOR, (1 - p_est[i]));  //（16）14
		amp_last[i] = (g[i] * amp_pr[i]);  //幅值
			if(i<5){

                DSP_X[i].r = 0*X[i].r;
                DSP_X[i].i = 0*X[i].i;
				//X[i].r *= 0;
				//X[i].i *= 0; 
		}
			else
			{
                DSP_X[i].r = g[i]*(X[i].r);
                DSP_X[i].i = g[i]*(X[i].i);
				//X[i].r  *= g[i] ;
				//X[i].i  *= g[i] ;
				}
		pr_snr[i] = MIN(pow(amp_last[i] / (noise_est[i] + 1e-30f), 2), 10000);  // 10000
		}
		else{
			X[i].r *= 0;
			X[i].i *= 0;
		}

	}
	
	return 0;
}
