clear 
clear all

figure
A1 = importdata('bubble.dat');
plot(A1.data(:,1),(A1.data(:,2)),'linewidth',2)
xlabel('Time(s)','FontName','Times New Roman','FontSize',15)
ylabel('Radius(m)','FontName','Times New Roman','FontSize',15)
set(gca,'LineWidth',1);
set(gca,'GridLineStyle','--','gridalpha',0.2,'FontName','Times New Roman','FontSize',15)
grid on

figure
A1 = importdata('bubble.dat');
plot(A1.data(:,1),(A1.data(:,3)),'linewidth',2)
xlabel('Time(s)','FontName','Times New Roman','FontSize',15)
ylabel('Migration(m)','FontName','Times New Roman','FontSize',15)
set(gca,'LineWidth',1);
set(gca,'GridLineStyle','--','gridalpha',0.2,'FontName','Times New Roman','FontSize',15)
grid on

figure
A2 = importdata('sensor.dat');
plot(A2.data(:,1),(A2.data(:,2)),'linewidth',2)
xlabel('Time(s)','FontName','Times New Roman','FontSize',15)
ylabel('Pressure(Pa)','FontName','Times New Roman','FontSize',15)
set(gca,'LineWidth',1);
set(gca,'GridLineStyle','--','gridalpha',0.2,'FontName','Times New Roman','FontSize',15)
grid on