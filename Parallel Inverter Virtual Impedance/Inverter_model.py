import numpy as np
pi = np.pi
Vdc = 700;
Lf, Rf, Cf,Rc1, Lc1,Rc2,Lc2 = 3.5E-3,0.1,50E-6,0.01,0.35E-3,0.05,0.8E-3;
Rv1,Rv2 = Rc1,Rc2;
Lv1,Lv2 = Lc1,Lc2;
def Inverter_Controller(t,x,RL,XL,m1,m2,n1,n2,Kpc,Kic,Kpv,Kiv):
    vdref = 400*np.sqrt(2/3)
    vqref = 0
    wf,w0 = 4*np.pi, 2*np.pi*50

    Pinv1,Qinv1,delta1,ild1,ilq1,vod1,voq1,iod1,ioq1,gammad1,gammaq1,zetad1,zetaq1,Pinv2,Qinv2,delta2,ild2,ilq2,vod2,voq2,iod2,ioq2,gammad2,gammaq2,zetad2,zetaq2 = x
    #Load equation
    ioD1, ioQ1 = -np.sin(delta1)*iod1 - np.cos(delta1)*ioq1, np.cos(delta1)*iod1 - np.sin(delta1)*ioq1
    ioD2, ioQ2 = -np.sin(delta2)*iod2 - np.cos(delta2)*ioq2, np.cos(delta2)*iod2 - np.sin(delta2)*ioq2
    vbD = (ioD1+ioD2)*RL - (ioQ1 + ioQ2)*XL
    vbQ = (ioD1 + ioD2)*XL + (ioQ1 + ioQ2)*RL
    vbd1,vbq1 = -np.sin(delta1)*vbD + np.cos(delta1)*vbQ, -np.cos(delta1)*vbD - np.sin(delta1)*vbQ
    vbd2,vbq2 = -np.sin(delta2)*vbD + np.cos(delta2)*vbQ, -np.cos(delta2)*vbD - np.sin(delta2)*vbQ
    #controller
    pinv1 = (3/2)*(vod1*iod1 + voq1*ioq1)
    qinv1 = (3/2)*(voq1*iod1 - vod1*ioq1)
    Pinv1_dot = - wf*Pinv1 + wf*pinv1 
    Qinv1_dot = - wf*Qinv1 + wf*qinv1 
    w1  = w0 - Pinv1 * m1
    delta1_dot = w1 - w1
    vodref1in = vdref - Qinv1*n1
    voqref1in = vqref;

    #Virtual impedance
    #vodref1in = u(8);
    #voqref1in = 0;
    #calulate vodref1 and voqref1
    vodref1 = vodref1in - iod1*Rv1 + ioq1*w1*Lv1;
    voqref1 = voqref1in - ioq1*Rv1 - iod1*w1*Lv1;

    ifdref1  = Kpv*(vodref1 - vod1) + Kiv*gammad1 + iod1 - w1*Cf*voq1
    ifqref1  = Kpv*(voqref1 - voq1) + Kiv*gammaq1 + ioq1 - w1*Cf*vod1
    gammad1_dot = vodref1 - vod1
    gammaq1_dot = voqref1 - voq1
    zetad1_dot = ifdref1 - ild1
    zetaq1_dot = ifqref1 - ilq1
    vfdref1 = Kpc*(ifdref1 - ild1) + Kic*zetad1 + vod1 - w1*Lf*ilq1
    vfqref1 = Kpc*(ifqref1 - ilq1) + Kic*zetaq1 + voq1 + w1*Lf*ild1

    vid1,viq1 = np.array([vfdref1,vfqref1])
    ild1_dot = (-Rf/Lf)*ild1 + w1*ilq1 + (vid1 - vod1)/Lf
    ilq1_dot = -w1*ild1 + (-Rf/Lf)*ilq1 + (viq1 - voq1)/Lf
    vod1_dot = w1*voq1 + (ild1 - iod1)/Cf
    voq1_dot = -w1*vod1 + (ilq1 - ioq1)/Cf
    iod1_dot = (-Rc1/Lc1)*iod1 + w1*ioq1 + (vod1 - vbd1)/Lc1
    ioq1_dot = -w1*iod1 + (-Rc1/Lc1)*ioq1 + (voq1 - vbq1)/Lc1


    #for inverter 2
    pinv2 = (3/2)*(vod2*iod2 + voq2*ioq2)
    qinv2 = (3/2)*(voq2*iod2 - vod2*ioq2)
    Pinv2_dot = - wf*Pinv2 + wf*pinv2
    Qinv2_dot = - wf*Qinv2 + wf*qinv2 
    w2 = w0 - Pinv2 * m2
    delta2_dot = w2 - w1
    vodref2in = vdref - Qinv2*n2
    voqref2in = vqref;
    #Virtual impedance
    #vodref1in = u(8);
    #voqref1in = 0;
    #calulate vodref1 and voqref1
    vodref2 = vodref2in - iod2*Rv2 + ioq2*w2*Lv2;
    voqref2 = voqref2in - ioq2*Rv2 - iod2*w2*Lv2;

    ifdref2  = Kpv*(vodref2 - vod2) + Kiv*gammad2 + iod2 - w2*Cf*voq2
    ifqref2  = Kpv*(voqref2 - voq2) + Kiv*gammaq2 + ioq2 - w2*Cf*vod2
    gammad2_dot = vodref2 - vod2
    gammaq2_dot = voqref2 - voq2
    zetad2_dot = ifdref2 - ild2
    zetaq2_dot = ifqref2 - ilq2
    vfdref2 = Kpc*(ifdref2 - ild2) + Kic*zetad2 + vod2 - w2*Lf*ilq2
    vfqref2 = Kpc*(ifqref2 - ilq2) + Kic*zetaq2 + voq2 + w2*Lf*ild2
    vid2,viq2 = np.array([vfdref2,vfqref2])

    ild2_dot = (-Rf/Lf)*ild2 + w2*ilq2 + (vid2 - vod2)/Lf
    ilq2_dot = -w2*ild2 + (-Rf/Lf)*ilq2 + (viq2 - voq2)/Lf
    vod2_dot = w2*voq2 + (ild2 - iod2)/Cf
    voq2_dot = -w2*vod2 + (ilq2 - ioq2)/Cf
    iod2_dot = (-Rc2/Lc2)*iod2 + w2*ioq2 + (vod2 - vbd2)/Lc2
    ioq2_dot = -w2*iod2 + (-Rc2/Lc2)*ioq2 + (voq2 - vbq2)/Lc2

    return np.array([Pinv1_dot,Qinv1_dot,delta1_dot,ild1_dot,ilq1_dot,vod1_dot,voq1_dot,iod1_dot,ioq1_dot,gammad1_dot,gammaq1_dot,zetad1_dot,zetaq1_dot,Pinv2_dot,Qinv2_dot,delta2_dot,ild2_dot,ilq2_dot,vod2_dot,voq2_dot,iod2_dot,ioq2_dot,gammad2_dot,gammaq2_dot,zetad2_dot,zetaq2_dot])

from sympy import symbols,Matrix,diff, sqrt,print_latex,zeros, sin,cos

# define symbols
#Inverter Symbols
ild1,ilq1,vod1,voq1,iod1,ioq1 = symbols(['i_ld1','i_lq1','v_od1','v_oq1','i_od1','i_oq1'])
ild2,ilq2,vod2,voq2,iod2,ioq2 = symbols(['i_ld2','i_lq2','v_od2','v_oq2','i_od2','i_oq2'])
wf,w0 = symbols(['omega_f','omega_0'])
vdref,vqref = symbols(['v_d_ref','v_q_ref'])
m1,m2,n1,n2 = symbols(['m_1','m_2','n_1','n_2'])
gammad1,gammad2,zetad1,zetad2,gammaq1,gammaq2,zetaq1,zetaq2 = symbols(['gamma_d1','gamma_d2','zeta_d1','zeta_d2','gamma_q1','gamma_q2','zeta_q1','zeta_q2'])
delta1,delta2,wf,Pinv1,Pinv2,Qinv1,Qinv2 = symbols(['delta1','delta2','w_f','P_inv1','P_inv2','Q_inv1','Q_inv2'])

def get_matrix(der,stat):
    mat = zeros(len(der),len(stat))
    for i in range(len(der)):
        for j in range(len(stat)):
            mat[i,j] = diff(der[i],stat[j])
    return mat

from sympy import symbols,Matrix,diff, sqrt,print_latex,zeros, sin,cos

def Inverter_Linear(RL,XL,m1,m2,n1,n2,Kpc,Kic,Kpv,Kiv):
        #Load equation
    ioD1, ioQ1 = -sin(delta1)*iod1 - cos(delta1)*ioq1, cos(delta1)*iod1 - sin(delta1)*ioq1
    ioD2, ioQ2 = -sin(delta2)*iod2 - cos(delta2)*ioq2, cos(delta2)*iod2 - sin(delta2)*ioq2
    vbD = (ioD1+ioD2)*RL - (ioQ1 + ioQ2)*XL
    vbQ = (ioD1 + ioD2)*XL + (ioQ1 + ioQ2)*RL
    vbd1,vbq1 = -sin(delta1)*vbD + cos(delta1)*vbQ, -cos(delta1)*vbD - sin(delta1)*vbQ
    vbd2,vbq2 = -sin(delta2)*vbD + cos(delta2)*vbQ, -cos(delta2)*vbD - sin(delta2)*vbQ

    #For inverter
    pinv1 = (3/4)*(vod1*iod1 + voq1*ioq1)
    qinv1 = (3/4)*(voq1*iod1 - vod1*ioq1)
    Pinv1_dot = - wf*Pinv1 + wf*pinv1 
    Qinv1_dot = - wf*Qinv1 + wf*qinv1 
    w1  = w0 - Pinv1 * m1
    delta1_dot = w1 - w1
    vodref1in = vdref - Qinv1*n1
    voqref1in = vqref;
    #Virtual impedance
    #vodref1in = u(8);
    #voqref1in = 0;
    #calulate vodref1 and voqref1
    vodref1 = vodref1in - iod1*Rv1 + ioq1*w1*Lv1;
    voqref1 = voqref1in - ioq1*Rv1 - iod1*w1*Lv1;

    ifdref1  = Kpv*(vodref1 - vod1) + Kiv*gammad1 + iod1 - w1*Cf*voq1
    ifqref1  = Kpv*(voqref1 - voq1) + Kiv*gammaq1 + ioq1 - w1*Cf*vod1
    gammad1_dot = vodref1 - vod1
    gammaq1_dot = voqref1 - voq1
    zetad1_dot = ifdref1 - ild1
    zetaq1_dot = ifqref1 - ilq1
    vfdref1 = Kpc*(ifdref1 - ild1) + Kic*zetad1 + vod1 - w1*Lf*ilq1
    vfqref1 = Kpc*(ifqref1 - ilq1) + Kic*zetaq1 + voq1 + w1*Lf*ild1

    vid1,viq1 = vfdref1,vfqref1
    ild1_dot = (-Rf/Lf)*ild1 + w1*ilq1 + (vid1 - vod1)/Lf
    ilq1_dot = -w1*ild1 + (-Rf/Lf)*ilq1 + (viq1 - voq1)/Lf
    vod1_dot = w1*voq1 + (ild1 - iod1)/Cf
    voq1_dot = -w1*vod1 + (ilq1 - ioq1)/Cf
    iod1_dot = (-Rc1/Lc1)*iod1 + w1*ioq1 + (vod1 - vbd1)/Lc1
    ioq1_dot = -w1*iod1 + (-Rc1/Lc1)*ioq1 + (voq1 - vbq1)/Lc1


    #for inverter 2
    pinv2 = (3/4)*(vod2*iod2 + voq2*ioq2)
    qinv2 = (3/4)*(voq2*iod2 - vod2*ioq2)
    Pinv2_dot = - wf*Pinv2 + wf*pinv2
    Qinv2_dot = - wf*Qinv2 + wf*qinv2 
    w2 = w0 - Pinv2 * m2
    delta2_dot = w2 - w1
    vodref2in = vdref - Qinv2*n2
    voqref2in = vqref;
    #Virtual impedance
    #vodref1in = u(8);
    #voqref1in = 0;
    #calulate vodref1 and voqref1
    vodref2 = vodref2in - iod2*Rv2 + ioq2*w2*Lv2;
    voqref2 = voqref2in - ioq2*Rv2 - iod2*w2*Lv2;

    ifdref2  = Kpv*(vodref2 - vod2) + Kiv*gammad2 + iod2 - w2*Cf*voq2
    ifqref2  = Kpv*(voqref2 - voq2) + Kiv*gammaq2 + ioq2 - w2*Cf*vod2
    gammad2_dot = vodref2 - vod2
    gammaq2_dot = voqref2 - voq2
    zetad2_dot = ifdref2 - ild2
    zetaq2_dot = ifqref2 - ilq2
    vfdref2 = Kpc*(ifdref2 - ild2) + Kic*zetad2 + vod2 - w2*Lf*ilq2
    vfqref2 = Kpc*(ifqref2 - ilq2) + Kic*zetaq2 + voq2 + w2*Lf*ild2
    vid2,viq2 = vfdref2,vfqref2

    ild2_dot = (-Rf/Lf)*ild2 + w2*ilq2 + (vid2 - vod2)/Lf
    ilq2_dot = -w2*ild2 + (-Rf/Lf)*ilq2 + (viq2 - voq2)/Lf
    vod2_dot = w2*voq2 + (ild2 - iod2)/Cf
    voq2_dot = -w2*vod2 + (ilq2 - ioq2)/Cf
    iod2_dot = (-Rc2/Lc2)*iod2 + w2*ioq2 + (vod2 - vbd2)/Lc2
    ioq2_dot = -w2*iod2 + (-Rc2/Lc2)*ioq2 + (voq2 - vbq2)/Lc2

    Xinv = Matrix([    Pinv1,     Qinv1,   delta1,     ild1,   ilq1,     vod1,    voq1,    iod1,   ioq1,   gammad1,     gammaq1,    zetad1,    zetaq1,    Pinv2,    Qinv2,    delta2,    ild2,     ilq2,    vod2,   voq2,    iod2,    ioq2,    gammad2,     gammaq2,   zetad2,    zetaq2])
    Xinv_dot = Matrix([Pinv1_dot,Qinv1_dot,delta1_dot,ild1_dot,ilq1_dot,vod1_dot,voq1_dot,iod1_dot,ioq1_dot,gammad1_dot,gammaq1_dot,zetad1_dot,zetaq1_dot,Pinv2_dot,Qinv2_dot,delta2_dot,ild2_dot,ilq2_dot,vod2_dot,voq2_dot,iod2_dot,ioq2_dot,gammad2_dot,gammaq2_dot,zetad2_dot,zetaq2_dot])
    Yinv = Matrix([Pinv1,Pinv2,Qinv1,Qinv2])
    Uinv = Matrix([vdref,vqref,w0])
    Ainv,Binv = get_matrix(Xinv_dot,Xinv),get_matrix(Xinv_dot,Uinv)
    #Cinv = get_matrix(Yinv,Xinv)
    return Ainv,Binv

def Linear_Inverter_Model(RL,XL,m1,m2,n1,n2,Kpc,Kic,Kpv,Kiv,X0):
    w0i,wfi = np.pi*2*50,np.pi*4
    Ainv,Binv= Inverter_Linear(RL,XL,m1,m2,n1,n2,Kpc,Kic,Kpv,Kiv)
    pinv1si,qinv1si,delta1si,ild1si,ilq1si,vod1si,voq1si,iod1si,ioq1si,gammad1si,gammaq1si,zetad1si,zetaq1si,pinv2si,qinv2si,delta2si,ild2si,ilq2si,vod2si,voq2si,iod2si,ioq2si,gammad2si,gammaq2si,zetad2si,zetaq2si = X0

    Amm = Ainv.subs([(iod1,iod1si),(iod2,iod2si),(ioq1,ioq1si),(ioq2,ioq2si),(vod1,vod1si),(voq1,voq1si),(vod2,vod2si),(voq2,voq2si)])
    Amm = Amm.subs([(w0,w0i),(wf,wfi),(delta1,delta1si),(delta2,delta2si),(Pinv1,pinv2si),(Pinv2,pinv2si)])
    Amm = Amm.subs([(ild1,ild1si),(ild2,ild2si),(ilq1,ilq1si),(ilq2,ilq2si)])

    Bmm = Binv.subs([(ild1,ild1si),(ild2,ild2si),(ilq1,ilq1si),(ilq2,ilq2si)])
    Bmm = Bmm.subs([(iod1,iod1si),(iod2,iod2si),(ioq1,ioq1si),(ioq2,ioq2si),(vod1,vod1si),(voq1,voq1si),(vod2,vod2si),(voq2,voq2si)])
    return np.array(Amm).astype(np.float64),np.array(Bmm).astype(np.float64)
