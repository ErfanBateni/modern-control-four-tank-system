# modern-control-four-tank-system
The objective of this project is to stabilize and control the liquid levels in the lower tanks of a four-tank system. The project covers the derivation of the nonlinear differential equations, linearization around an equilibrium point, and the design and comparative analysis of various modern control techniques.

# Modern Control: Four-Tank System Analysis & Control

This repository contains the Phase 2 project for the **Modern Control** course at Sharif University of Technology. The project focuses on the modeling, analysis, and controller design for a **Quadruple-Tank Process**, a well-known benchmark for non-linear, multi-input multi-output (MIMO) systems exhibiting complex dynamics such as cross-coupling and parameter variations.

## 📌 Project Overview
The objective of this project is to stabilize and control the liquid levels in the lower tanks of a four-tank system. The project covers the derivation of the nonlinear differential equations, linearization around an equilibrium point, and the design and comparative analysis of various modern control techniques.

Key topics covered include:
*   **System Modeling:** Linearization of the nonlinear MIMO dynamics into a state-space representation and diagonalization using the Jordan canonical form.
*   **State-Feedback Control:** Pole-placement design aiming for specific transient performance (settling time < 60s, overshoot < 40%).
*   **Precompensator Design:** Implementation of static and dynamic (integral action) precompensators for reference tracking and constant disturbance rejection.
*   **State Observers:** Design and comparison of Full-Order and Reduced-Order Luenberger observers under noisy measurement conditions.
*   **Optimal Control (LQR):** Tuning Q and R matrices for balanced energy consumption and transient performance.
*   **Robustness Analysis:** Evaluating the performance of observers and controllers (Pole-placement vs. LQR) under a 20% variation in the system's state matrix ($A$).

## ⚙️ Key Findings & Results
*   **Linear vs. Nonlinear Simulation:** State-feedback controllers designed for the linearized model perform exceptionally well near the operating point (e.g., $t_s \approx 46s$, overshoot $\approx 21\%$). However, when applied to the full nonlinear model, performance degrades significantly (overshoot > 65%, steady-state error) due to saturation limits and nonlinearity, highlighting the need for gain scheduling or nonlinear Model Predictive Control (MPC).
*   **Disturbance Rejection:** Dynamic precompensators successfully reject constant disturbances ($y_{ss} = 0$), whereas static precompensators leave a steady-state offset.
*   **Observer Performance:** The Full-Order observer provides smoother estimations and better noise attenuation. The Reduced-Order observer is computationally lighter but highly sensitive to measurement noise on the unmeasured states.
*   **Robustness:** Under a 20% system variation, the LQR controller maintained robust, non-oscillatory stability, whereas the pole-placement controller induced complex conjugate poles leading to oscillations. Redesigning the controllers based on the updated plant model successfully restored nominal performance.

## 🛠️ Tools & Technologies
*   **MATLAB & Simulink:** State-space modeling, control design (`place`, `lqr`, `lsim`), and numerical simulations.
*   **Modern Control Theory:** MIMO systems, Luenberger Observers, Linear Quadratic Regulator (LQR), Robustness Analysis.

## 📂 Repository Structure
*   `/src`: Contains all MATLAB scripts (`Q1.m` to `Q9.m`) corresponding to each section of the analysis (Linearization, Observers, LQR, Robustness, etc.).
*   `Report_Phase2_400100792.pdf`: The comprehensive technical report (in Persian) detailing mathematical derivations, control strategies, and plotted simulation results.

## 👥 Author
*   **Erfan Bateni (400100792)**
