from ECMOGUI import Ui_MainWindow
from headloss import HeadLoss as HL
import sys
# import matplotlib
# matplotlib.use("Qt5Agg")  # 声明使用QT5
import numpy as np
# import matplotlib.pyplot as plt
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui, QtWidgets
# from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
# from matplotlib.figure import Figure
# 全局变量
g=9.8
cmTom=0.01
cmTomm=10
mmTom=0.001
inchtocm=2.54
FrTomm=1/3
Ltom3=0.001
mintos=60
mmHgtoPa=133.322
PatommHg=1/133.322
# 定义窗口类，继承Ui_MainWindow
class myMainWindow(QtWidgets.QMainWindow,Ui_MainWindow):
# 初始化，继承父类所有参数
    def __init__(self):
        super(myMainWindow,self).__init__()
        self.setupUi(self)
# 菜单按钮的相应事件
        self.actionExit.triggered.connect(self.on_click_exit)
        self.actionExit.setShortcut('Ctrl+Q')
        self.actionExit.setStatusTip('Exit application')
# 程序中计算和退出按钮的响应事件
        self.pushButton.clicked.connect(self.on_click_compute)
        self.pushButton_2.clicked.connect(self.on_click_clean)
# 点击菜单栏按钮的响应
    def on_click_exit(self):
        sys.exit(app.exec_())
# 清除前一次计算结果
    def on_click_clean(self):
        self.textBrowser.clear()
        self.textBrowser_2.clear()
        self.textBrowser_3.clear()
        self.textBrowser_4.clear()
        self.textBrowser_5.clear()
        self.textBrowser_6.clear()
# 点击计算时的响应
    def on_click_compute(self):
# 读取所有输入数据的信息
        try:
            D_p = float(self.lineEdit.text())
            D_vc = float(self.lineEdit_2.text())
            D_ac = float(self.lineEdit_3.text())
            L_vc = float(self.lineEdit_4.text())
            L_ac = float(self.lineEdit_5.text())

            Q_in = float(self.lineEdit_6.text())
            H_a= float(self.lineEdit_7.text())
            P_outv= float(self.lineEdit_8.text())
            P_outa= float(self.lineEdit_9.text())

            rho= float(self.lineEdit_10.text())
            mu= float(self.lineEdit_11.text())

            theta_b= float(self.lineEdit_12.text())
            L_vp= float(self.lineEdit_13.text())
            L_ap= float(self.lineEdit_14.text())
            theta_vp= float(self.lineEdit_15.text())
            theta_ap= float(self.lineEdit_16.text())
            L_vt= float(self.lineEdit_17.text())
            L_at= float(self.lineEdit_18.text())
            self.textBrowser_4.insertPlainText("数据读取成功...\n")
            # 调用单位转换
            No1=Preprocess(D_p,D_vc,D_ac,L_vc,L_ac,Q_in,H_a,P_outv,P_outa,rho,mu,theta_b,L_vp,L_ap,theta_vp,theta_ap,L_vt,L_at)
            D_p,D_vc,D_ac,L_vc,L_ac,Q_in,H_a,P_outv,P_outa,rho,mu,theta_b,L_vp,L_ap,theta_vp,theta_ap,L_vt,L_at=No1.Unit_transform()
            self.textBrowser_4.insertPlainText("单位转换成功...\n")
        except:
            self.textBrowser_4.insertPlainText("数据读取失败，请检查后重试\n")


# 调用计算模块
        try:
            No2=Solve(D_p,D_vc,D_ac,L_vc,L_ac,Q_in,H_a,P_outv,P_outa,rho,mu,theta_b,L_vp,L_ap,theta_vp,theta_ap,L_vt,L_at)
            beta_v=No2.compute()
            self.textBrowser_4.insertPlainText("求解成功...\n")
# 调用后处理模块
            No3=Postprocess(beta_v,Q_in)
            Q_v,Q_a=No3.Dwdata()
# 将结果显示回GUI
            self.textBrowser.insertPlainText(str(round(Q_v,4))+"\n")
            self.textBrowser_2.insertPlainText(str(round(Q_a,4))+"\n")
            self.textBrowser_3.insertPlainText(str(round(beta_v,4))+"\n")
            self.textBrowser_5.insertPlainText(str(round(beta_v*Q_in/A_vc,4))+"\n")
            self.textBrowser_6.insertPlainText(str(round((1-beta_v)*Q_in/A_ac,4))+"\n")
        except:
            # 调用后处理模块
            beta_v=1.0
            No3=Postprocess(beta_v,Q_in)
            Q_v,Q_a=No3.Dwdata()
            self.textBrowser_4.insertPlainText("设置错误或流量超过1，未求得流量或单向状态，请重试\n")
            self.textBrowser.insertPlainText(str(Q_v)+"\n")
            self.textBrowser_2.insertPlainText(str(Q_a)+"\n")
            self.textBrowser_3.insertPlainText(str(beta_v)+"\n")
            self.textBrowser_5.insertPlainText(str(beta_v*Q_in/A_vc)+"\n")
            self.textBrowser_6.insertPlainText(str((1-beta_v)*Q_in/A_ac)+"\n")

# 前处理模块，包括数据的单位处理等
class Preprocess(object):
    def __init__(self,D_p,D_vc,D_ac,L_vc,L_ac,Q_in,H_a,P_outv,P_outa,rho,mu,theta_b,L_vp,L_ap,theta_vp,theta_ap,L_vt,L_at):
        self.D_p=D_p
        self.D_vc=D_vc
        self.D_ac=D_ac
        self.L_vc=L_vc
        self.L_ac=L_ac
        self.Q_in=Q_in
        self.H_a=H_a
        self.P_outv=P_outv
        self.P_outa=P_outa
        self.rho=rho
        self.mu=mu
        self.theta_b=theta_b
        self.L_vp=L_vp
        self.L_ap=L_ap
        self.theta_vp=theta_vp
        self.theta_ap=theta_ap
        self.L_vt=L_vt
        self.L_at=L_at
# 单位的变换
    def Unit_transform(self):
        self.D_p=self.D_p*inchtocm*cmTom
        self.D_vc=(self.D_vc-2)*FrTomm*mmTom
        self.D_ac=(self.D_ac-2)*FrTomm*mmTom
        self.L_vc=self.L_vc*cmTom
        self.L_ac=self.L_ac*cmTom
        self.Q_in=self.Q_in*Ltom3/mintos
        self.H_a=self.H_a*cmTom
        self.P_outv=self.P_outv*mmHgtoPa
        self.P_outa=self.P_outa*mmHgtoPa
        self.rho=self.rho
        self.mu=self.mu
        self.L_vp=self.L_vp*cmTom
        self.L_ap=self.L_ap*cmTom
        self.theta_vp=self.theta_vp
        self.theta_ap=self.theta_ap
        self.L_vt=self.L_vt*cmTom
        self.L_at=self.L_at*cmTom
        return self.D_p,self.D_vc,self.D_ac,self.L_vc,self.L_ac,self.Q_in,self.H_a,self.P_outv,self.P_outa,self.rho,self.mu,self.theta_b,self.L_vp,self.L_ap,self.theta_vp,self.theta_ap,self.L_vt,self.L_at
# 求解计算模块
class Solve(object):
    def __init__(self,D_p,D_vc,D_ac,L_vc,L_ac,Q_in,H_a,P_outv,P_outa,rho,mu,theta_b,L_vp,L_ap,theta_vp,theta_ap,L_vt,L_at):
        self.D_p=D_p
        self.D_vc=D_vc
        self.D_ac=D_ac
        self.L_vc=L_vc
        self.L_ac=L_ac
        self.Q_in=Q_in
        self.H_a=H_a
        self.P_outv=P_outv
        self.P_outa=P_outa
        self.rho=rho
        self.mu=mu
        self.theta_b=theta_b
        self.L_vp=L_vp
        self.L_ap=L_ap
        self.theta_vp=theta_vp
        self.theta_ap=theta_ap
        self.L_vt=L_vt
        self.L_at=L_at
    def compute(self):
        theta_vt=np.rad2deg(np.arctan((self.D_p-self.D_vc)/(2*self.L_vt)))
        theta_at=np.rad2deg(np.arctan((self.D_p-self.D_ac)/(2*self.L_at)))
        global A_p,A_vc,A_ac
        nu=self.mu/self.rho
        A_p=self.D_p**2*np.pi/4
        A_vc=self.D_vc**2*np.pi/4
        A_ac=self.D_ac**2*np.pi/4
        beta0=0
        beta=0.5
        beta1=1
        residual=1
        U_in=self.Q_in/A_p
        n=1
    # 二分法求解
        while np.abs(residual)>0.0001 and n<100:
            n+=1
            #连续方程
            U_vp=beta*U_in
            U_ap=U_in-U_vp
            U_vc=U_vp*A_p/A_vc
            U_ac=(U_in-U_vp)*A_p/A_ac
            #弯管损失
            H_Labe=HL.ElbowLoss((360-self.theta_b)/2,U_ap)
            H_Lvbe=HL.ElbowLoss((360-self.theta_b)/2,U_vp)

            H_Lve=HL.ElbowLoss(180-self.theta_vp,U_vp)
            H_Lae=HL.ElbowLoss(180-self.theta_ap,U_ap)
            #渐缩管损失
            H_Lvt=HL.TaperLoss(theta_vt,self.D_p,self.D_vc,U_vc)
            H_Lat=HL.TaperLoss(theta_at,self.D_p,self.D_ac,U_ac)
            #导管沿程损失
            H_fvp=HL.FlowLoss(self.L_vp,self.D_p,U_vp,nu)
            H_fap=HL.FlowLoss(self.L_ap,self.D_p,U_ap,nu)
            #插管沿程损失
            H_fvc=HL.FlowLoss(self.L_vc,self.D_vc,U_vc,nu)
            H_fac=HL.FlowLoss(self.L_ac,self.D_ac,U_ac,nu)
            #分岔损失系数计算
            U=np.array([U_in,-U_vp,-U_ap])
            A=np.array([A_p,A_p,A_p])
            theta=np.array([-np.pi/2,3*np.pi/4,np.pi/4])
            k=HL.JunctionLossCoefficient(U,A,theta)
            H01v=0.5*k[0]*U_in**2/g
            H01a=0.5*k[1]*U_in**2/g
            #动量方程
            P_vp=self.P_outv+0.5*self.rho*(U_vc**2-U_vp**2)+self.rho*g*(H_fvc+H_fvp+H_Lvt+H_Lve+H_Lvbe)
            P_ap=self.P_outa+0.5*self.rho*(U_ac**2-U_ap**2)+self.rho*g*(self.H_a+H_fac+H_fap+H_Lae+H_Labe+H_Lat)
            P_inv=P_vp+0.5*self.rho*(U_vp**2-U_in**2)+self.rho*g*H01v
            P_ina=P_ap+0.5*self.rho*(U_ap**2-U_in**2)+self.rho*g*H01a

            residual=P_inv-P_ina
            if residual>0:
                jie=beta
                beta=(beta0+beta)/2
                beta1=jie
            else:
                jie=beta
                beta=(beta+beta1)/2
                beta0=jie 
        if n>=99:
            beta=1.0
        return beta
 #后处理模块
class Postprocess(object):
    def __init__(self,beta_v,Q_in):
        self.beta_v=beta_v
        self.Q_in=Q_in
    def Dwdata(self):
        self.Q_in=self.Q_in*mintos/Ltom3
        Q_v=self.beta_v*self.Q_in
        Q_a=(1-self.beta_v)*self.Q_in
        return Q_v,Q_a

# 主函数
if __name__=='__main__':
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app=QApplication(sys.argv)
    ui=myMainWindow()
    ui.show()
    sys.exit(app.exec_())