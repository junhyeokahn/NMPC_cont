
figure()
subplot(2,1,1)
plot(t,state_nom(:,7:8)*(180/pi),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('Attitude [deg]');
legend('Roll','Pitch');

subplot(2,1,2)
plot(t,ctrl_dyn_nom(:,2:4)*(180/pi),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('Angular rate [deg/s]');
hl = legend('$\dot{\phi}$','$\dot{\theta}$','$\dot{\psi}$');
set(hl,'interpreter','latex');

figure()
subplot(2,1,1)
plot(t, state_nom(:,1:3),'linewidth',2);
grid on
xlabel('Time [s]'); legend('x','y','z');

subplot(2,1,2)
plot(t, state_nom(:,4:6),'linewidth',2);
grid on
xlabel('Time [s]'); legend('v_x','v_y','v_z');

figure()
subplot(2,1,1)
plot(t, accel_nom,'linewidth',2);
grid on
xlabel('Time [s]'); legend('a_x','a_y','a_z');      

subplot(2,1,2)
plot(t, ctrl_dyn_nom(:,1),'linewidth',2);
grid on
xlabel('Time [s]'); legend('Thrust derivative');