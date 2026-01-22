%Input
L2 = 40; %cm
L3 = 30; %cm
Omega2 = 30; %rad/s
Omega3 = 20; %rad/s
t = 0:0.001:2;
IniTeta2 = pi/4;
IniTeta3 = pi/6;
NumData = size(t,2);

q = zeros(9,1);
q_all = zeros(9,NumData);
qdot_all = zeros(9,NumData);
qdot2_all = zeros(9,NumData);

for j = 1:NumData
    for i = 1:10 %jumlah iterasi
        C = [q(1); q(2); q(3);...
            q(4)-(L2*cos(q(6))/2);...
            q(5)-(L2*sin(q(6))/2);...
            q(4)+(L2*cos(q(6))/2)-q(7)+(L3*cos(q(9))/2);...
            q(5)+(L2*sin(q(6))/2)-q(8)+(L3*sin(q(9))/2);...
            q(6)-IniTeta2-Omega2*t(j);...
            q(9)-IniTeta3-Omega3*t(j)];

        %Jacobian Matrix
        Cq = [1 0 0 0 0 0 0 0 0;
              0 1 0 0 0 0 0 0 0;
              0 0 1 0 0 0 0 0 0;
              0 0 0 1 0 (L2*sin(q(6))/2) 0 0 0;
              0 0 0 0 1 -(L2*cos(q(6))/2) 0 0 0;
              0 0 0 1 0 -(L2*sin(q(6))/2) -1 0 -(L3*sin(q(9))/2);
              0 0 0 0 1 (L2*cos(q(6))/2) 0 -1 (L3*cos(q(9))/2);
              0 0 0 0 0 1 0 0 0;
              0 0 0 0 0 0 0 0 1];

        q_delta = - inv(Cq)*C;
        q = q + q_delta;

        %velocity Analysis
        Ct = [0;0;0;0;0;0;0;-Omega2;-Omega3;];
        qdot = -inv(Cq)*Ct;

        %Acceleration Analysis
        Cq_qdot_q = [0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 Omega2*L2*cos(q(6))/2 0 0 0;
            0 0 0 0 0 Omega2*L2*sin(q(6))/2 0 0 0;
            0 0 0 0 0 -Omega2*L2*cos(q(6))/2 0 0 -Omega3*L3*cos(q(9))/2;
            0 0 0 0 0 -Omega2*L2*sin(q(6))/2 0 0 -Omega3*L3*sin(q(9))/2;
            0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0];

        Qd = -Cq_qdot_q*qdot;
        qdot2 = inv(Cq)*Qd;

    end

    q_all(:,j) = q;
    qdot_all(:,j) = qdot;
    qdot2_all(:,j) = qdot2;

end
figure(1)
plot(t,q_all(7,:),t,q_all(8,:));


figure(2)
plot(q_all(7,:),q_all(8,:));

figure(3)
plot(q_all(4,:),q_all(5,:));

figure(3)
plot(q_all(4,:),q_all(5,:));

figure(4)
plot(q_all(3,:),q_all(6,:));


