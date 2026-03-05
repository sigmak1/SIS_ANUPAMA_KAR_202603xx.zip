# SIS_ANUPAMA_KAR_202603xx.zip

**Class:** AE 544, Analytical Dynamics  
**Professor:** Dr. Hao Peng  
**Semester:** Spring 2026  
**Student Name:** Anupama Kar   
**Assignment:** Programming Project 01   
**School:** Embry-Riddle Aeronautical University

--------------  
Usage of AI: Used AI to understand and learn how to set everything up on GitHub, to help in markdown language,  MathJax formatting for equations, and to learn how to create .gif from MATLAB, to polish my base code, to define some of the functions, and to add some fancy items like indicator plots. Also used it to research and read about different integrators to enhance my 
knowledge base. However, I did verify the results using other references as stated below, in the Reference section. 

## Table of Contents
1. [How-To-Guide](#1-how-to-guide)
2. [Singularity for asymmetric Euler angle sets](#2-singularity-for-asymmetric-euler-angle-sets)
   - [2.1. Background and set-up of the EOM](#21-background-and-set-up-of-the-eom)
   - [2.2. Results](#22-results)
   - [2.3. Analysis](#23-analysis)
3. [Ambiguity of Euler parameters (quaternions)](#3-ambiguity-of-euler-parameters-quaternions)
   - [3.1. Background and set-up of the EOM](#31-background-and-set-up-of-the-eom)
   - [3.2. Results](#32-results)
   - [3.3. Analysis](#33-analysis)

4. [Classical Rodrigues Parameters (CRP)](#4-classical-rodrigues-parameters-crp)
   - [4.1. Background and set-up of the EOM](#41-background-and-set-up-of-the-eom)
   - [4.2. Results](#42-results)
   - [4.3. Analysis](#43-analysis)
5. [Comparison of different numerical integrators](#5-comparison-of-different-numerical-integrators)
   - [5.1. Overview of ODE Solvers Used](#51-overview-of-ode-solvers-used)
   - [5.2. Result](#52-result)
   - [5.3. Analysis](#53-analysis)
6. [3D animation explanation](#6-3d-animation-explanation)
7. [References](#7-references)

## List of Figures
Figure 2.1: Successive yaw, pitch, and roll rotations  
Figure 2.2: 3-2-1 Euler Angles vs Time (ode45)  
Figure 3.1: Quaternion Components   
Figure 3.2: Quaternion Norm  
Figure 4.1: CRP Norm and Principal Rotation Angle  
Figure 5.1: Integrator Step Size vs Time vs Pitch  
Figure 5.2: Singularity syntax error encountered during simulation  
Figure 6.1: Aircraft Attitude (.gif)

## List of Abbreviations
EOM: Equation of Motion  
DCM: Direction Cosine Matrix  
CRP: Classical Rodrigues Parameters  
ODE: Ordinary Differential Equation  
MATLAB: Matrix Laboratory (software)  
GIF: Graphics Interchange Format

## 1. How-To-Guide
Follow the below step by step guide to navigate through the folder:

1. Log on to your Github account. 
2. Find the zip folder: https://github.com/sigmak9/SIS_ANUPAMA_KAR_202603xx.zip 
3. Download and extract all.
4. Open *README.md*
5. Open *main.m*. This is written on MATLAB R2025b version.   
**Note:** The functions, and the code was written using a very useful Youtube video by Dr. Shane Ross, from Virginia Tech [1], and then Gen AI was used to clean up. 
6. Launch MATLAB and copy the *main.m* script to it and run the script.
7. Let the animation run for 60 seconds. Do not click on any other plot while the animation is running as that can interfere with the results.  

Basically, the animation code it written using Gen AI and has not been verified but the way it works is, if you click on a different figure while the animation is running, it will lose that figure and will activate the animation on the other figure.   

8. Once the simulation hits t=60s, you should have given 6 figures (including the animation dialogue box).
9. The animation will be saved on the same directory/ path as a *Aircraft_Attitude.gif* file. 
10. The remaining figures/ plots are as follows:  
   a. Figure 1: 3-2-1 Euler Angles vs Time (ode45)   
   b. Figure 2: Integrator Step Size vs Time and Pitch ($\theta$)  
   c. Figure 3: Quaternion Components (Look for Sign Flip)  
   d. Figure 4: Quaternion Norm (Should Remain 1)  
   e: Figure 5: CRP Norm and Principal Rotation Angle 

The same single MATLAB script will generate all the required plots and results for the asks. The Euler Angles, Euler Parameters, Classical Rodeiguez Parameters (CRPs), three integrator comparisons, and some indicator plots. Basically **Figures 4** is a good indicator plot to make sure the computation is done right.

## 2.	Singularity for asymmetric Euler angle sets
### 2.1. Background and set-up of the EOM  
Euler angles describe the orientation of a rigid body with respect to a fixed coordinate system using three successive rotations. Let $(\phi, \theta, \psi)$ be the three Euler angles (often called roll, pitch, and yaw).

![Figure 2.1: Successive yaw, pitch, and roll rotations](Successive%20yaw,%20pitch,%20and%20rotations.png)

**Figure 2.1:** Successive yaw, pitch, and roll rotations. This figure illustrates the sequence of rotations about the principal axes, corresponding to the Euler angles $(\psi, \theta, \phi)$.

The transformation from the body frame to the inertial frame can be represented by the rotation matrix $R$:

$$
R = M_3(\psi)\,M_2(\theta)\,M_1(\phi)
$$

$$
M_1(\phi) =
\begin{bmatrix}
1 & 0 & 0 \\
0 & \cos\phi & \sin\phi \\
0 & \sin\phi & \cos\phi
\end{bmatrix}
$$

$$
M_2(\theta) =
\begin{bmatrix}
\cos\theta & 0 & \sin\theta \\
0 & 1 & 0 \\
\sin\theta & 0 & \cos\theta
\end{bmatrix}
$$

$$
M_3(\psi) =
\begin{bmatrix}
\cos\psi & \sin\psi & 0 \\
\sin\psi & \cos\psi & 0 \\
0 & 0 & 1
\end{bmatrix}
$$


The body angular velocity vector in general form is:

$$
\omega_B = \begin{bmatrix} \omega_x \\ \omega_y \\ \omega_z \end{bmatrix}
$$

In the MATLAB script, the specific angular velocity vector is defined as a function of time $t$:

```matlab
w_fun = @(t) [0.2*sin(0.02*t);
				0.15;
				0.2*cos(0.02*t)];
```

Designate the Direction Cosine Matrix (DCM) relating Euler angle rates to body angular velocity as matrix $B$:

$$
B = \begin{bmatrix}
1 & 0 & -\sin\theta \\
0 & \cos\phi & \sin\phi \cos\theta \\
0 & -\sin\phi & \cos\phi \cos\theta
\end{bmatrix}
$$

Assign the Euler angle rates as:

$$
\dot{y} = \begin{bmatrix} \dot{\psi} \\ \dot{\theta} \\ \dot{\phi} \end{bmatrix}
$$  

**Equation (3.57)** [2] shows the Euler Angles kinematic equation of motion as reproduced below:
<div style="border:2px solid black; padding:10px; display:inline-block;">
The Euler Angles kinematic equation of motion is:

$$
\dot{y} = B^{-1} \omega_B
$$
</div>

where 

$$
B^{-1} = \frac{1}{\cos \theta} 
\begin{bmatrix}
0 & \sin \phi & \cos \phi \\
0 & \cos \theta \cos \phi & -\cos \theta \sin \phi \\
\cos \theta & \sin \theta \sin \phi & \sin \theta \cos \phi
\end{bmatrix}
$$

**Note: Here B in the code is equal to $B^{-1}$ in the derivation.**
 
In the *main.m* script, it is coded as a function 

```matlab
function dydt = EulerODE(t,y,w_fun)
% y = [psi theta phi]'
psi = y(1);
theta = y(2);
phi = y(3);

w = w_fun(t);

s2 = sin(theta);
c2 = cos(theta);
s3 = sin(phi);
c3 = cos(phi);

B = 1/c2*[0 s3 c3;
          0 c2*c3 -c2*s3;
          c2 s2*s3 s2*c3];

dydt = B*w;

end
```
Looking at the inverse B matrix, the **singularity occurs when $\cos \theta = 0$, i.e., $\theta = \pi/2$ or 90°**. At this point, the factor $1/\cos\theta$ becomes **infinite**, making the matrix undefined. 

Physically, this corresponds to **gimbal lock**, where the Euler angles lose one degree of freedom: the rotation axes align such that you cannot uniquely determine all three angular rates from the Euler angle rates. Near $\theta = 90^\circ$, small changes in angular velocity cause large changes in the Euler rates, making control or simulation unstable.

### 2.2. Results

![Figure 2.2: 3-2-1 Euler Angles vs Time (ode45)](3-2-1%20Euler%20Angles%20ve%20Time%20(ode45).png)

**Figure 2.2:** 3-2-1 Euler Angles vs Time (ode45). This figure shows the time evolution of the Euler angles using the ode45 integrator.

### 2.3. Analysis  
**Euler singularity is geometric** (depends on the middle angle).

The pitch angle θ steadily increases and approaches 90° as designed in the angular velocity profile (**Figure 2.1**). Near t ≈ 25–30 s, θ gets very close to 90°, signaling the approach to the Euler singularity. At this point, the other two angles, yaw ψ and roll φ, begin to show rapid, erratic changes or large jumps. When θ → 90°, the rotation matrix becomes nearly singular, and ψ and φ are no longer uniquely defined. In the plot, this appears as step-ups or step-downs in yaw and roll — sudden changes in their values while the physical motion remains smooth. This is a classic visual cue for Euler angle singularities: the angles “jump” because the coordinate mapping cannot represent rotations uniquely near θ = ±90°.

Singularities can be spotted in the plot by observing when the pitch angle θ approaches ±90° (dashed reference lines are included in the plot) and by noting where ψ or φ suddenly increase or decrease sharply. Often, adaptive integrators like `ode45` or `ode15s` reduce their step size near the singularity, which can also appear as denser points along the curve. Fixed-step integrators such as `ode4` do not adapt and may produce larger numerical errors in yaw and roll near the singularity.

## 3.	Ambiguity of Euler parameters (quaternions)  
### 3.1. Background and set-up of the EOM 
To curb the singularity issue in Euler angles, four Euler parameters are used, mainly $\beta_0$, $\beta_1$, $\beta_2$, and $\beta_3$ (or $q_0$, $q_1$, $q_2$, and $q_3$).

<div style="border:2px solid black; padding:10px; display:inline-block;">
The Euler Parameters kinematic equation of motion is:

$$
\dot{\mathbf{q}} = \frac{1}{2} \, \Omega \, \mathbf{q}, \quad  
\Omega =
\begin{bmatrix}
0 & \omega_x & \omega_y & \omega_z \\
-\omega_x & 0 & -\omega_z & \omega_y \\
-\omega_y & \omega_z & 0 & -\omega_x \\
-\omega_z & -\omega_y & \omega_x & 0
\end{bmatrix}
$$
</div>  

Here, $\dot{\mathbf{q}}$ = $\dot{\boldsymbol{\beta}}$, $\mathbf{q}$ = $\boldsymbol{\omega}_B$, and $B(\boldsymbol{\beta})$ = $\Omega$.

The Junkins and Schaub textbook gives a thorough explanation and derivation of the quarternions. **Equation (3.104)** [2] is shown below:

$$
\dot{\boldsymbol{\beta}} =
\frac{1}{2}
\begin{bmatrix}
0 & -\omega_1 & -\omega_2 & -\omega_3 \\
\omega_1 & 0 & \omega_3 & -\omega_2 \\
\omega_2 & -\omega_3 & 0 & \omega_1 \\
\omega_3 & \omega_2 & -\omega_1 & 0
\end{bmatrix}
\begin{bmatrix}
\beta_0 \\ \beta_1 \\ \beta_2 \\ \beta_3
\end{bmatrix}
$$

where 
The quaternion column vector is:

$$
\boldsymbol{\beta} =
\begin{bmatrix}
\beta_0 \\
\beta_1 \\
\beta_2 \\
\beta_3
\end{bmatrix}
$$
or 
$$
\mathbf{q} =
\begin{bmatrix}
q_0 \\
q_1 \\
q_2 \\
q_3
\end{bmatrix}
$$  

This is defined as a function on the *main.m* file as 

````matlab
function dqdt = QuatODE(t,q,w_fun)
% q = [b0 b1 b2 b3]'
w = w_fun(t);

Omega = [ 0    -w(1) -w(2) -w(3);
          w(1)  0     w(3) -w(2);
          w(2) -w(3)  0     w(1);
          w(3)  w(2) -w(1)  0];

dqdt = 0.5 * Omega * q;

end
````
Quaternions have a **sign ambiguity** because a quaternion and its negative represent the same physical rotation. That is, $\mathbf{q} = [q_0, q_1, q_2, q_3]^T$ and $-\mathbf{q} = [-q_0, -q_1, -q_2, -q_3]^T$ produce identical rotations in 3D space. This can cause apparent “jumps” in quaternion values when integrating or interpolating rotations if the sign flips, even though the actual orientation does not change. Care must be taken in numerical algorithms to maintain sign consistency to avoid discontinuities.

### 3.2. Results  

![Figure 3.1: Quaternion Components](Quaternion%20Components.png)

**Figure 3.1:** Quaternion components. This figure shows the individual components of the quaternion over time.

![Figure 3.2: Quaternion Norm](Quaternion%20Norm.png)

**Figure 3.2:** Quaternion norm. This figure shows the norm of the quaternion, which should remain close to 1 for a valid rotation.

### 3.3. Analysis  

The quaternion-based integration produces several key plots: the individual components b₀, b₁, b₂, b₃, the quaternion norm ‖q‖, and the dot product between consecutive quaternions to detect ambiguity.

1. **Quaternion Components vs Time:**  
   The plots of b₀, b₁, b₂, b₃ show smooth evolution even as the spacecraft approaches large rotations or Euler singularities. Occasionally, one component may appear to "flip" sign. This happens because of the **quaternion double-cover property**: q and -q represent the same physical rotation.  

   Visually, a flip is seen when a decreasing curve becomes increasing abruptly.   

   The first noticeable flip occurs around **t ≈ 12.6 s** in **Figure 3.1**. Importantly, these flips are **coordinate artifacts**, not physical discontinuities.

2. **Quaternion Norm:**  
   The plot of ‖q‖ = √(b₀² + b₁² + b₂² + b₃²) remains extremely close to 1 throughout the simulation (**Figure 3.2**). This confirms that the quaternions are properly normalized and immune to the singularities that affect Euler angles or CRP. Normalization is explicitly enforced in the code via the `normalizeQuat` function. This is a good indicator plot to proof read all math, calculation and computation.

## 4.	Classical Rodrigues Parameters (CRP)  
### 4.1. Background and set-up of the EOM  
<div style="border:2px solid black; padding:10px; display:inline-block;">
The standard CRP kinematic equation is:

$$
\dot{\mathbf{r}} = \frac{1}{2} \left( \mathbf{I}_3 + [\mathbf{r}]_\times + \mathbf{r} \mathbf{r}^T \right) \boldsymbol{\omega}_B
$$
</div>

where

$$
[\mathbf{r}]_\times =
\begin{bmatrix}
0 & -r_3 & r_2 \\
r_3 & 0 & -r_1 \\
-r_2 & r_1 & 0
\end{bmatrix}.
$$

This matches the form given on **Equation (3.131)** [2]. 

On the MATLAB file, this CRP function is defined as below

````matlab
function drdt = CRPODE(t,r,w_fun)

w = w_fun(t);

r1 = r(1);
r2 = r(2);
r3 = r(3);

r_tilde = [  0   -r3   r2;
             r3    0   -r1;
            -r2   r1    0];

drdt = 0.5 * ( eye(3) + r_tilde + r*r' ) * w;

end
````

CRP was developed to mitigate the ambiguity issue arising from the quaternions. 
Classical Rodrigues Parameters (CRP) are defined from a quaternion 
$\mathbf{q} = [q_0, \mathbf{q}_v]^T$, where $\mathbf{q}_v = [q_1, q_2, q_3]^T$, as:

$$
\mathbf{r} = \frac{\mathbf{q}_v}{q_0} =
\frac{1}{q_0}
\begin{bmatrix}
q_1 \\ q_2 \\ q_3
\end{bmatrix}, \quad q_0 \neq 0 \ (\text{i.e., rotation angle } \theta \neq 180^\circ)
$$

The mapping from a rotation to CRP is **unique**, except at $\theta = 180^\circ$.

While this solves the ambiguity issue, it reintroduces the singularity issue at $\theta = 180^\circ$.  

### 4.2. Results  

![Figure 4.2: CRP Norm and Principal Rotation Angle](CRP%20Norm%20and%20Principal%20Rotation%20Angle.png)

**Figure 4.1:** CRP Norm and Principal Rotation Angle. This figure shows the relationship between the CRP norm and the principal rotation angle, highlighting the singularity as the angle approaches $180^\circ$.

### 4.3. Analysis  
The combined plot of the CRP norm and principal rotation angle as shown in **Figure 4.1** shows that the classical Rodrigues parameters become unbounded as the principal rotation angle approaches $180^\circ$. While the quaternion-based principal angle evolves smoothly through $180^\circ$, the CRP norm spikes because

$$
|\mathbf{r}| = \tan\left(\frac{\Phi}{2}\right),
$$

which tends to infinity as $\Phi \to 180^\circ$. This demonstrates that the CRP singularity depends on the total rotation angle $\Phi$ rather than any individual Euler angle component.

As the rotation approaches $180^\circ$, the CRP state grows extremely large, causing the numerical integrator step size to shrink and eventually fail because the state becomes unbounded. Consequently, no meaningful CRP solution is produced beyond this point. The aircraft attitude itself remains continuous and well-defined, but the CRP parameterization becomes invalid, confirming that this behavior is a geometric singularity of the representation rather than a physical discontinuity in the motion.

**CRP singularity is global** (depends on total rotation). 

For this kinematics, the total rotation hits $180^\circ$ at time **t=~12.59 sec**. 

## 5.	Comparison of different numerical integrators    
### 5.1. Overview of ODE Solvers Used

In this project, three different numerical integration methods are employed to solve the attitude kinematics: `ode45`, `ode15s`, and `ode4`. The built-in MATLAB solvers such as `ode45` and `ode15s` are adaptive, variable-step methods that adjust step size based on the local behavior of the solution to satisfy specified error tolerances.  

 `ode45` implements a medium-order explicit Runge–Kutta (Dormand–Prince) algorithm and is generally effective for most nonstiff ordinary differential equations, making it a good first choice for typical dynamics problems.  
 
  `ode15s` is a variable-step, variable-order solver based on numerical differentiation formulas, designed to handle stiff systems or problems that are inefficient or difficult for `ode45` to integrate.  
  
   In contrast, `ode4` is a simple fixed-step fourth-order Runge–Kutta integrator with a user-specified step size; it does not automatically adapt to rapid changes in the solution because it lacks error control. These differences in algorithm and adaptivity influence how each solver behaves as the attitude kinematics approach singular configurations. [3]

### 5.2. Result
![Figure 5.1: Integrator Step Size vs Time vs Pitch](Integrator%20Step%20Size%20vs%20Time%20vs%20Pitch.png)

**Figure 5.1:** Integrator step size vs time vs pitch. This figure illustrates how the integrator step size varies with time and pitch angle during the simulation.  

![Singularity syntax error](Singularity%20syntax%20error.png)

**Figure 5.2:** Singularity syntax error encountered during simulation. This illustrates the numerical breakdown as the CRP singularity is reached.  

### 5.3. Analysis

The step-size plot (**Figure 5.1**) shows distinct behavior for `ode45`, `ode15s`, and the fixed-step `ode4` method as the Euler singularity (θ → ±90°) or the CRP singularity (Φ → 180°) is approached.

As the singularity is approached, the kinematic equations become ill-conditioned and the derivatives grow rapidly. For `ode45` (an explicit adaptive Runge–Kutta 4/5 method), the solver automatically reduces its step size to maintain the specified error tolerances. This appears in the plot as a sharp decrease in step size near the singularity. Although `ode45` can approach the singularity smoothly, its step size may become extremely small, increasing computational cost and potentially causing termination if the required step becomes too small.

`ode15s` (a variable-order stiff solver based on numerical differentiation formulas) is designed for stiff systems. Near the singularity, the rapid growth in derivatives introduces stiffness-like behavior. As a result, `ode15s` typically maintains stability better than `ode45` and may require fewer total steps, but it can still struggle when the state itself becomes unbounded (as in CRP at Φ = 180°). Since the singularity is geometric rather than purely stiffness-driven, `ode15s` does not eliminate the blow-up; it only handles the rapid dynamics more robustly.

In contrast, `ode4` uses a fixed step size and does not adapt to the rapidly increasing derivatives. As the singularity is approached, its fixed step may be too large to accurately resolve the steep gradients, leading to large numerical errors or instability. The step-size plot for `ode4` remains constant, but the solution accuracy degrades significantly near the singularity.

Overall, `ode45` provides the best balance of accuracy and efficiency for non-stiff attitude kinematics, as it adapts its step size appropriately. `ode15s` may offer improved stability in rapidly changing regions but does not resolve the fundamental singularity. The fixed-step `ode4` method is the least reliable near singularities because it cannot adjust its step size to maintain accuracy.

The solver failure warning from `ode45` at **t ≈ 12.6 s** as shown in **Figure 5.2** coincides with the principal rotation angle approaching 180° (**Figure 4.1**). As the CRP norm diverges, the state derivatives grow without bound, forcing the adaptive solver to reduce its step size below machine precision. This numerical breakdown confirms that the observed behavior is due to the geometric singularity of the CRP parameterization rather than a physical instability in the spacecraft motion. This warning is strong evidence that the CRP singularity has been reached.

## 6. 3D animation explanation  
Once MATLAB runs successfully, it will generate a .gif file and save it in the designated folder as *Aircraft_Attitude.gif*. The legends of the .gif are same as the legends of the Euler Angles plot (figure 1.2.). Since the .gif updates itiratively, putting a legend that would stay fixed was rendering bad graphics.  


![Figure 6.1: Aircraft Attitude](Aircraft_Attitude.gif)

**Figure 6.1:** Aircraft attitude. This figure shows the aircraft attitude animation as a GIF.

It is to be noted that the animation is a rendering of the computation and in no way represents the actual physical motion of an aircraft. So, singularity in Euler angles or CRPs results in tumbling or rapid spinning. Quarternion ambiguity does not show any symptoms. Only singularity affects the actual computer attitude.  

I was usign AI and trying to tinker around with putting models into the animation as shown in [1]. But I could only put a 2D picture so we will deal with this for now. Also, the Yaw, Pitch, and Roll vectors do not line up with the aircrafts nose, right wing, and upward as it should. 
## 7. References
[1] Ross Dynamics Lab, “Euler Angle Simulation with MATLAB – Integrating the Rotational Kinematic Differential Equations,” YouTube video, 30 Jul. 2021. https://www.youtube.com/watch?v=vwn_JT0SDXQ. 

[2] Schaub, H., and Junkins, J. L., Analytical Mechanics of Space Systems, 2nd ed., American Institute of Aeronautics and Astronautics, Reston, VA, 2009.  

[3] MathWorks, “Choose an ODE Solver,” MATLAB & Simulink Documentation, The MathWorks, Inc., Natick, MA, accessed Mar. 2, 2026. Available: https://www.mathworks.com/help/matlab/math/choose-an-ode-solver.html