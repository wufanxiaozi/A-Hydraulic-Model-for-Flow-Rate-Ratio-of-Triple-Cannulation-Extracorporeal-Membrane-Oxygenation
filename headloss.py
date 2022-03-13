import numpy as np
g=9.8        
Re1=2040
Re2=3200
class HeadLoss(object):
    #局部损失函数
    def FlowLoss(L,d,U,nu):
        Re=U*d/nu
        #lam = 0.3164/Re**(1/4) 湍流
        #lam = 64/Re 层流
        #层流、湍流均适用
        op=operator()
        a,b=op.TurnAround(Re1,Re2)
        if Re<Re1:
            lam=64/Re
        elif Re<=Re2:
            lam=a*Re**b
        else:
            lam=0.3164/Re**(1/4)
        H = L/d*U**2/(2*g)*lam
        return H


    def ElbowLoss(theta,U):
        H=(0.945*np.sin(np.deg2rad((180-theta)/2))**2+2.047*np.sin(np.deg2rad((180-theta)/2))**4)*0.5*U**2/g
        return H

    def TaperLoss(theta,D1,D2,U):
        if D1<D2:
            error("不是渐缩管")
        if theta<=22.5:
            H=0.8*np.sin(np.deg2rad(theta))*(1-(D2/D1)**2)/(D2/D1)**4*U**2/(2*g)
        elif theta<=90:
            H=0.5*(1-(D2/D1)**2)*np.sqrt(np.sin(np.deg2rad(theta)))/(D2/D1)**4*U**2/(2*g)
        else:
            error("该角度不在范围内")
        return H

    def JunctionLossCoefficient(U,A,theta):
        U=U.ravel()
        A=A.ravel()
        op=operator()
        theta = op.wrapToPi([theta.ravel()])
        theta = theta.ravel()
        Q = U*A
        Ci = Q < 0.0
        Si = Q >= 0.0
        Qtot = np.sum(Q[Si])
        FlowRatio = -Q[Ci]/Qtot
        PseudoColAngle = np.mean(theta[~Si])
        PseudoSupAngle = np.arctan2(np.sum(np.sin(theta[Si])*Q[Si]), np.sum((np.cos(theta[Si]))*Q[Si]))
        if np.abs(PseudoSupAngle-PseudoColAngle) < np.pi/2:
            PseudoColAngle = PseudoColAngle + np.pi
        theta = op.wrapToPi(theta - PseudoColAngle)

# Calculate the pseudosupplier angle
        pseudodirection = np.sign(np.mean(np.sin(theta[Si])*Q[Si])) 
        if pseudodirection < 0:  
            theta = -theta
        PseudoSupAngle = np.arctan2(np.sum(np.sin(np.abs(theta[Si]))*Q[Si]), np.sum(np.cos(np.abs(theta[Si]))*Q[Si]))  # always positive
# Calculate effective pseudosupplier area
        etransferfactor = (0.8*(np.pi-PseudoSupAngle)*np.sign(theta[Ci])-0.2)*(1-FlowRatio)
        TotPseudoArea = Qtot / ((1 - etransferfactor)*np.sum(U[Si]*Q[Si])/Qtot)  # Denominator is pseudosupplier velocity  

# Calculate area ratios and relative angles
        AreaRatio = TotPseudoArea / A[Ci];
        phi = op.wrapTo2Pi(PseudoSupAngle - theta[Ci])

# Calculate the C loss coefficients
        C = np.zeros(U.shape[0])
        C[Ci] = (1-np.exp(-FlowRatio/0.02))*(1-(1.0/(AreaRatio*FlowRatio))*np.cos(0.75*(np.pi-phi)))

# Calculate the K loss coefficients (if required)
        if U.shape[0] <= 3:
            if np.sum(Ci) == 1:
                Ucom = U[Ci]
            else:
                Ucom = U[Si]
            K = (U[Ci]**2/Ucom**2)*(2*C[Ci] + U[Si]**2/U[Ci]**2-1)
        return K

class operator(object):
    def wrapToPi(self,theta):
        thetawrap=np.remainder(theta, 2*np.pi)
        mask = np.abs(thetawrap)>=np.pi
        thetawrap[mask] -= 2*np.pi * np.sign(thetawrap[mask])
        return thetawrap
    def wrapTo2Pi(self,theta):
        thetawrap=np.remainder(theta, 2*np.pi)
        mask = np.abs(thetawrap)>=2*np.pi
        thetawrap[mask] -= 2*np.pi * np.sign(thetawrap[mask])
        return thetawrap
    def TurnAround(self,Re1,Re2):
#   转捩区域模型
#   给定转捩开始雷诺数和结束雷诺数，
#   返回 $\lambda$ =a*Re^b的系数（双对数下是线性增加）
        b=np.log10(64*Re2**(0.25)/(Re1*0.3164))/np.log10(Re1/Re2)
        a=64*Re1**(-1-b)
        return a,b