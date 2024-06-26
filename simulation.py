from sympy import symbols, cos, sin, pi
from sympy.matrices import Matrix

import matplotlib.pyplot as plt
import numpy as np

a1=1
a2=1
a3=5
a4=5
a5=0.2


#----------------------------------------------------completed-------------------------------

def DH_Params_H_Matrix(i,j,q1=0,q2=0,q3=0,q4=0,q5=0,q6=0):
    DH_params=np.array([
        [q1,np.pi/2,a2,a1],
        [q2,0,a3,0],
        [q3+np.pi/2,np.pi/2,0,0],
        [q4,-np.pi/2,0,a4],
        [q5,np.pi/2,0,0],
        [q6,0,0,a5]
    ])
    H_matrix= Matrix([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])
    for k in range(i,j):
        H2_matrix=Matrix([
            [np.cos(DH_params[k,0]),-np.sin(DH_params[k,0])*np.cos(DH_params[k,1]),np.sin(DH_params[k,0])*np.sin(DH_params[k,1]),DH_params[k,2]*np.cos(DH_params[k,0])],
            [np.sin(DH_params[k,0]),np.cos(DH_params[k,1])*np.cos(DH_params[k,0]),-np.cos(DH_params[k,0])*np.sin(DH_params[k,1]),DH_params[k,2]*np.sin(DH_params[k,0])],
            [0,np.sin(DH_params[k,1]),np.cos(DH_params[k,1]),DH_params[k,3]],
            [0,0,0,1]
        ])
        H_matrix=np.dot(H_matrix,H2_matrix)
    return np.matrix(H_matrix)


#----------------------------------------------------completed-------------------------------

def forward_kinematics(q1,q2,q3,q4,q5,q6):
    X=[0,0]
    Y=[0,0]
    Z=[0,a1]
    for k in range(1,7):
        H_matrix=DH_Params_H_Matrix(0,k,q1,q2,q3,q4,q5,q6)
        X.append(H_matrix[0,3])
        Y.append(H_matrix[1,3])
        Z.append(H_matrix[2,3])
    X = np.reshape(X, (1, 8))
    Y = np.reshape(Y, (1, 8))
    Z = np.reshape(Z, (1, 8))
    return   X,Y,Z
    


#----------------------------------------------------completed-------------------------------

def cosinelaw(a,b,c):
    return np.arccos((a*a + b*b - c*c) / (2*a*b))
     


#----------------------------------------------------completed-------------------------------

   

def Gripper_Coordinates(finalX,finalY,finalZ,R0_6):
    Xgripper=finalX-a5*R0_6[0,2]
    Ygripper=finalY-a5*R0_6[1,2]
    Zgripper=finalZ-a5*R0_6[2,2]
    return Xgripper,Ygripper,Zgripper
    


#----------------------------------------------------completed-------------------------------


def inverse_kinemtics_first_3(Xgripper,Ygripper,Zgripper):
    #print(Xgripper,Ygripper,Zgripper)
    q1=np.arctan2(Ygripper,Xgripper)
    r1=np.sqrt(Xgripper**2+Ygripper**2)
    r2=np.sqrt((r1-a2)**2+(Zgripper-a1)**2)
    phi1=np.arctan((Zgripper-a1)/(r1-a2))
    phi2=cosinelaw(a3,r2,a4)
    q2=phi1+phi2
    phi3=cosinelaw(a3,a4,r2)
    q3=phi3-np.pi
    return q1,q2,q3
    


#----------------------------------------------------completed-------------------------------


def inverse_kinemtics_last_3(R3_6):
    
    q5=np.arccos(R3_6[2,2])
    q6=np.arcsin(R3_6[2,1]/np.sin(q5))
    q4=np.arcsin(R3_6[1,2]/np.sin(q5))
    return q4,q5,q6
    


#----------------------------------------------------completed-------------------------------


def get_Joint_Variables(finalX,finalY,finalZ,roll,pitch,yaw):
    R0_6_initial=np.matrix([
        [0,0,1],
        [0,-1,0],
        [1,0,0]
    ])
    R_6_rotation=np.matrix([
        [np.cos(pitch)*np.cos(roll),np.sin(pitch)*np.sin(yaw)-np.sin(roll)*np.cos(pitch)*np.cos(yaw),np.sin(pitch)*np.cos(yaw)+np.sin(roll)*np.sin(yaw)*np.cos(pitch)],
        [np.sin(roll),np.cos(roll)*np.cos(yaw),-np.sin(yaw)*np.cos(roll)],
        [-np.sin(pitch)*np.cos(roll),np.sin(pitch)*np.sin(roll)*np.cos(yaw)+np.sin(yaw)*np.cos(pitch),-np.sin(pitch)*np.sin(roll)*np.sin(yaw)+np.cos(pitch)*np.cos(yaw)]
    ])
    
    R0_6=np.dot(R0_6_initial,R_6_rotation)
    #print(R0_6)
    
    XGripper,YGripper,ZGripper=Gripper_Coordinates(finalX,finalY,finalZ,R0_6)
    q1,q2,q3=inverse_kinemtics_first_3(XGripper,YGripper,ZGripper)
    #print(q1,q2,q3)
    
    H0_3=DH_Params_H_Matrix(0,3,q1,q2,q3)
    R0_3=H0_3[0:3,0:3]
    
    
    Inverse_R0_3=np.linalg.inv(np.matrix(R0_3,dtype=np.float64))
    
    R3_6=np.dot(Inverse_R0_3,R0_6)
    #print(np.matrix(R3_6))
    #print(type(R3_6))
    q4,q5,q6=inverse_kinemtics_last_3(R3_6)
    #print(q1,q2,q3,q4,q5,q6)
    return q1,q2,q3,q4,q5,q6
    


#----------------------------------------------------completed-------------------------------

def create_plot():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('x axis')
    ax.set_ylabel('y axis')
    ax.set_zlabel('z axis')
    ax.set_autoscale_on(False)
    # fig.canvas.draw()
    # plt.show()
    return fig, ax

def update_plot( finalX,finalY,finalZ,X1, Y1, Z1,  fig, ax):
    GripperEndX = np.reshape(finalX, (1, 8))
    GripperEndY = np.reshape(finalY, (1, 8))
    GripperEndZ = np.reshape(finalZ, (1, 8))
    ax.cla()
    ax.plot_wireframe(GripperEndX, GripperEndY, GripperEndZ)
    ax.plot_wireframe(X1, Y1, Z1, color = 'g')
    # print('in update ', X1, Y1, Z1)
    plt.draw()
    ax.set_xlabel('x axis')
    ax.set_ylabel('y axis')
    ax.set_zlabel('z axis')
    ax.set_autoscale_on(False)
    fig.canvas.draw()
    fig.canvas.flush_events()
    plt.pause(0.5)

def simulate_word():
    fig, ax = create_plot()

    # corner points of word
    X0 = [4,4,4,4,4]
    Y0 = [10,10,9,8,8]
    Z0 = [1,3,2,3,1]
    X, Y, Z, X1, Y1, Z1 = [], [], [], [], [], []
    for i in range(5):
        #finalX,finalY,finalZ ARE THE COORDINATES OF THE GRIPPER END THAT IS AT NEW LOCATION ON EVERY JOINT OF LETTER
        finalX,finalY,finalZ = X0[i], Y0[i], Z0[i]
        X1 = X0[:i+1]; Y1 = Y0[:i+1]; Z1 = Z0[:i+1]
        # X1,Y1,Z1 ARE FOR PLOTTING WIREFRAME OF THE LETTER THAT HAS TO BE PRINTED
        X1 = np.reshape(X1, (1, i+1))
        Y1 = np.reshape(Y1, (1, i+1))
        Z1 = np.reshape(Z1, (1, i+1))
        # print(X1, Y1, Z1)
        roll, pitch, yaw = 0, 0, 0

        #NEEDED TO GET THE ROTATION OF EVERY SERVO AND TO FURTHER CALCULATE THE POSITION OF EVERY SERVO IN SPACE
        q1,q2,q3,q4,q5,q6=get_Joint_Variables(finalX,finalY,finalZ,roll,pitch,yaw)

        #ALL COORDINATES OF THE SERVOS JOINT ARRAY OF 8 POINTS TO DRAW WIREFRAME OF ROBOTIC ARM
        X,Y,Z=forward_kinematics(q1,q2,q3,q4,q5,q6)
        update_plot(X, Y, Z, X1, Y1, Z1,  fig, ax)

    # print(X1, Y1, Z1)
    # update_plot(X, Y, Z, X1, Y1, X1, fig, ax)
    plt.show()


#----------------------------------------MAIN PART---------------------------------
#finalX,finalY,finalZ=5,1,1
#roll,pitch,yaw=0,0,0
#q1,q2,q3,q4,q5,q6=get_Joint_Variables(finalX,finalY,finalZ,roll,pitch,yaw)
#X,Y,Z=forward_kinematics(q1,q2,q3,q4,q5,q6)
#print(X,Y,Z,sep='\n')
def main():
    simulate_word()
if __name__ == "__main__":
    main()