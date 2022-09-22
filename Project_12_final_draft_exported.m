classdef Project_12_final_draft_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        ObliqueIncidenceLabel           matlab.ui.control.Label
        WaveParametersLabel             matlab.ui.control.Label
        Ei0EditFieldLabel               matlab.ui.control.Label
        Ei0EditField                    matlab.ui.control.NumericEditField
        EditFieldLabel                  matlab.ui.control.Label
        EditField                       matlab.ui.control.NumericEditField
        iEditFieldLabel                 matlab.ui.control.Label
        iEditField                      matlab.ui.control.NumericEditField
        r1EditFieldLabel                matlab.ui.control.Label
        r1EditField                     matlab.ui.control.NumericEditField
        r1EditField_2Label              matlab.ui.control.Label
        r1EditField_2                   matlab.ui.control.NumericEditField
        r2EditFieldLabel                matlab.ui.control.Label
        r2EditField                     matlab.ui.control.NumericEditField
        r2EditField_2Label              matlab.ui.control.Label
        r2EditField_2                   matlab.ui.control.NumericEditField
        PerpendicularPolarizationLabel  matlab.ui.control.Label
        PlotE1EiErButton                matlab.ui.control.Button
        PlotE2EtButton                  matlab.ui.control.Button
        PlotH1HiHrButton                matlab.ui.control.Button
        PlotH2HtButton                  matlab.ui.control.Button
        BrewsterAngleBEditFieldLabel    matlab.ui.control.Label
        BrewsterAngleBEditField         matlab.ui.control.NumericEditField
        ParallelPolarizationLabel       matlab.ui.control.Label
        PlotE1EiErButton_2              matlab.ui.control.Button
        PlotE2EtButton_2                matlab.ui.control.Button
        PlotH1HiHrButton_2              matlab.ui.control.Button
        PlotH2HtButton_2                matlab.ui.control.Button
        CalculateBrewsterAngleButton    matlab.ui.control.Button
        CalculateBrewsterAngleButton_2  matlab.ui.control.Button
        BrewsterAngleBEditField_2Label  matlab.ui.control.Label
        BrewsterAngleBEditField_2       matlab.ui.control.NumericEditField
        degreeLabel                     matlab.ui.control.Label
        degreeLabel_2                   matlab.ui.control.Label
        degreeLabel_3                   matlab.ui.control.Label
        Label                           matlab.ui.control.Label
        MediumParametersLabel           matlab.ui.control.Label
        PlotE1E2Button                  matlab.ui.control.Button
        PlotE1E2Button_2                matlab.ui.control.Button
        PlotH1H2Button                  matlab.ui.control.Button
        PlotH1H2Button_2                matlab.ui.control.Button
    end

    
    
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.Ei0EditField.Value = 2;
            app.EditField.Value = 2*pi*10e8;
            app.iEditField.Value = 60;
            app.r1EditField.Value = 5;
            app.r1EditField_2.Value = 3;
            app.r2EditField.Value = 6;
            app.r2EditField_2.Value = 4;
        end

        % Button pushed function: CalculateBrewsterAngleButton
        function CalculateBrewsterAngleButtonPushed(app, event)
            
            epsilon1 = app.r1EditField.Value;
            mu1 = app.r1EditField_2.Value;
            epsilon2 = app.r2EditField.Value;
            mu2 = app.r2EditField_2.Value;
            eta1 = 120*pi*sqrt(mu1/epsilon1);
            %eta1
            eta2 = 120*pi*sqrt(mu2/epsilon2);
            %eta2
            theta_br = asind(sqrt((1 - (eta1/eta2)^2)/(1 - (mu1/mu2)^2)));
            app.BrewsterAngleBEditField.Value = real(theta_br);
            im_theta_br = imag(theta_br);
            if mu1==mu2
                app.Label.Text='Value of Brewster Angle does not exist in case of Perpendicular Polarization'
            elseif im_theta_br~=0
                app.Label.Text='Value does not exist as Brewster angle is imaginary'
            else
                app.Label.Text='';
            end
        end

        % Button pushed function: PlotE1EiErButton
        function PlotE1EiErButtonPushed(app, event)
            Ei0 = app.Ei0EditField.Value;
            w = app.EditField.Value;
            theta_i = app.iEditField.Value;
            epsilon1 = app.r1EditField.Value;
            mu1 = app.r1EditField_2.Value;
            epsilon2 = app.r2EditField.Value;
            mu2 = app.r2EditField_2.Value;
            
            eta1 = 120*pi*sqrt(mu1/epsilon1);
            %eta1
            eta2 = 120*pi*sqrt(mu2/epsilon2);
            %eta2
            e0 = 8.8541878128*10^-12;
            m0 = 4*pi*10^-7;
            beta1 = w*sqrt(mu1*epsilon1*e0*m0);
            %beta1
            beta2 = w*sqrt(mu2*epsilon2*e0*m0);
            %beta2
            theta_r = theta_i;
            %theta_r
            theta_t = asind(sind(theta_i)*beta1/beta2);
            %theta_t
            gamma1 = ((eta2*cosd(theta_i)) - (eta1*cosd(theta_t)))/((eta2*cosd(theta_i)) + (eta1*cosd(theta_t)));
            %gamma1
            x = linspace(-0.5,0,101);
            z = linspace(-0.5,0,101);
            y = linspace(-0.5,0,101);
            %defining E field Variable
            Eix = zeros(101,101);
            Eiy = zeros(101,101);
            Eiz = zeros(101,101);
            Erx = zeros(101,101);
            Ery = zeros(101,101);
            Erz = zeros(101,101);
            %setting values
            for n = 1:1:101
                for m = 1:1:101
                    Eiy(n, m) = Ei0.*exp(-1i*(beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                    Ery(n, m) = Ei0.*gamma1.*exp(-1i*(-beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                end
            end
            subplot(2,2,1);
            surf(x,z,real(Eiy));
            title('Perpendicular Polarization Electric Field');
            subtitle('Ei');
            xlabel('X-axis');
            ylabel('Z-axis');
            zlabel('Y-axis');
            subplot(2,2,2);
            surf(x,z,real(Ery));
            
            subtitle('Er');
            xlabel('X-axis');
            ylabel('Z-axis');
            zlabel('Y-axis');
            subplot(2,2,[3,4]);
            surf(x,z,real(Ery + Eiy));
            
            subtitle('Total E1(Ei + Er)');
            xlabel('X-axis');
            ylabel('Z-axis');
            zlabel('Y-axis');
            
        end

        % Button pushed function: PlotE2EtButton
        function PlotE2EtButtonPushed(app, event)
            Ei0 = app.Ei0EditField.Value;
            w = app.EditField.Value;
            theta_i = app.iEditField.Value;
            epsilon1 = app.r1EditField.Value;
            mu1 = app.r1EditField_2.Value;
            epsilon2 = app.r2EditField.Value;
            mu2 = app.r2EditField_2.Value;
            
            eta1 = 120*pi*sqrt(mu1/epsilon1);
            %eta1
            eta2 = 120*pi*sqrt(mu2/epsilon2);
            %eta2
            ep = 8.8541878128*10^-12;
            m = 4*pi*10^-7;
            beta1 = w*sqrt(mu1*epsilon1*ep*m);
            %beta1
            beta2 = w*sqrt(mu2*epsilon2*ep*m);
            %beta2
            
            theta_t = asind(sind(theta_i)*beta1/beta2);
            %theta_t
            
            tau1 = (2*eta2*cosd(theta_i))/((eta2*cosd(theta_i)) + (eta1*cosd(theta_t)));
            %tau1
            
            %Plot for E2
            x = linspace(0,0.5,101);
            z = linspace(0,0.5,101);
            y = linspace(0,0.5,101);
            %defining E field Variable
            Etx = zeros(101,101);
            Ety = zeros(101,101);
            Etz = zeros(101,101);
            
            %setting values
            for n = 1:1:101
                for m = 1:1:101
                    Ety(n, m) = tau1.*Ei0.*exp(-1i*(beta2.*cosd(theta_t).*x(n) + beta2.*sind(theta_t).*z(m)));
                    
                end
            end
            surf(x,z,real(Ety));
            title('Perpendicular Polarization Electric Field');
            subtitle('E2(Et)');
            xlabel('X-axis');
            ylabel('Z-axis');
            zlabel('Y-axis');
        end

        % Button pushed function: PlotH1HiHrButton
        function PlotH1HiHrButtonPushed(app, event)
            Ei0 = app.Ei0EditField.Value;
            w = app.EditField.Value;
            theta_i = app.iEditField.Value;
            epsilon1 = app.r1EditField.Value;
            mu1 = app.r1EditField_2.Value;
            epsilon2 = app.r2EditField.Value;
            mu2 = app.r2EditField_2.Value;
            
            eta1 = 120*pi*sqrt(mu1/epsilon1);
            %eta1
            eta2 = 120*pi*sqrt(mu2/epsilon2);
            %eta2
            e0 = 8.8541878128*10^-12;
            m0 = 4*pi*10^-7;
            beta1 = w*sqrt(mu1*epsilon1*e0*m0);
            %beta1
            beta2 = w*sqrt(mu2*epsilon2*e0*m0);
            %beta2
            theta_r = theta_i;
            %theta_r
            theta_t = asind(sind(theta_i)*beta1/beta2);
            %theta_t
            gamma1 = ((eta2*cosd(theta_i)) - (eta1*cosd(theta_t)))/((eta2*cosd(theta_i)) + (eta1*cosd(theta_t)));
            %gamma1
            x = linspace(-0.5,0,101);
            z = linspace(-0.5,0,101);
            y = linspace(-0.5,0,101);
            Hix = zeros(101, 101);
            Hiy = zeros(101, 101);
            Hiz = zeros(101, 101);
            
            Hrx = zeros(101, 101);
            Hry = zeros(101, 101);
            Hrz = zeros(101, 101);
            for n = 1:1:101
                for m = 1:1:101
                    Hix(n, m) = -sind(theta_i).*(Ei0/eta1).*exp(-1i*(beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                    Hiz(n, m) = cosd(theta_i).*(Ei0/eta1).*exp(-1i*(beta1.*cosd(theta_i).*x(n) +  beta1.*sind(theta_i).*z(m)));
                    Hrx(n, m) = -sind(theta_i).*(gamma1.*Ei0/eta1).*exp(-1i*(-beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                    Hrz(n, m) = -cosd(theta_i).*(gamma1.*Ei0/eta1).*exp(-1i*(-beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                end
            end
            subplot(2,2,1);
            surf(x, z, real(Hix + Hiz));
            title('Perpendicular Polarization Magnetic Field');
            subtitle('Hi');
            xlabel('X-axis');
            ylabel('Z-axis');
            zlabel('Y-axis');
            subplot(2,2,2);
            surf(x, z, real(Hrx + Hrz));
            
            subtitle('Hr');
            xlabel('X-axis');
            ylabel('Z-axis');
            zlabel('Y-axis');
            subplot(2,2,[3,4]);
            surf(x, z, real(Hix + Hrx + Hiz + Hrz));
            
            subtitle('Total H1(Hi + Hr)');
            xlabel("X-axis")
            ylabel("Z-axis")
            zlabel("Y-axis")
        end

        % Button pushed function: PlotH2HtButton
        function PlotH2HtButtonPushed(app, event)
            Ei0 = app.Ei0EditField.Value;
            w = app.EditField.Value;
            theta_i = app.iEditField.Value;
            epsilon1 = app.r1EditField.Value;
            mu1 = app.r1EditField_2.Value;
            epsilon2 = app.r2EditField.Value;
            mu2 = app.r2EditField_2.Value;
            
            eta1 = 120*pi*sqrt(mu1/epsilon1);
            %eta1
            eta2 = 120*pi*sqrt(mu2/epsilon2);
            %eta2
            ep = 8.8541878128*10^-12;
            m = 4*pi*10^-7;
            beta1 = w*sqrt(mu1*epsilon1*ep*m);
            %beta1
            beta2 = w*sqrt(mu2*epsilon2*ep*m);
            %beta2
            
            theta_t = asind(sind(theta_i)*beta1/beta2);
            %theta_t
            
            tau1 = (2*eta2*cosd(theta_i))/((eta2*cosd(theta_i)) + (eta1*cosd(theta_t)));
            %tau1
            
            %Plot for E2
            x = linspace(0,0.5,101);
            z = linspace(0,0.5,101);
            y = linspace(0,0.5,101);
            %defining E field Variable
            Htx = zeros(101,101);
            Hty = zeros(101,101);
            Htz = zeros(101,101);
            
            %setting values
            for n = 1:1:101
                for m = 1:1:101
                    Htx(n, m) = -sind(theta_t).*(tau1.*Ei0/eta2).*exp(-1i*(beta2.*cosd(theta_t).*x(n) + beta2.*sind(theta_t).*z(m)));
                    Htz(n, m) = cosd(theta_t).*(tau1.*Ei0/eta2).*exp(-1i*(beta2.*cosd(theta_t).*x(n) +  beta2.*sind(theta_t).*z(m)));
                end
            end
            surf(x, z, real(Htx + Htz));
            title('Perpendicular Polarization Magnetic Field');
            subtitle('Total H2(Ht)');
            xlabel('X-axis');
            ylabel('Z-axis');
            zlabel('Y-axis');
        end

        % Button pushed function: CalculateBrewsterAngleButton_2
        function CalculateBrewsterAngleButton_2Pushed(app, event)
            epsilon1 = app.r1EditField.Value;
            mu1 = app.r1EditField_2.Value;
            epsilon2 = app.r2EditField.Value;
            mu2 = app.r2EditField_2.Value;
            eta1 = 120*pi*sqrt(mu1/epsilon1);
            %eta1
            eta2 = 120*pi*sqrt(mu2/epsilon2);
            %eta2
            theta_br = asind(sqrt((1 - (eta2/eta1)^2)/(1 - (epsilon1/epsilon2)^2)));
            app.BrewsterAngleBEditField_2.Value = real(theta_br);
            im_theta_br = imag(theta_br);
            if epsilon1==epsilon2
                app.Label.Text='Value of Brewster Angle does not exist in case of Parallel Polarization'
            elseif im_theta_br~=0
                app.Label.Text='Value does not exist as Brewster angle is imaginary'
            else
                app.Label.Text='';
            end
        end

        % Button pushed function: PlotE1EiErButton_2
        function PlotE1EiErButton_2Pushed(app, event)
            Ei0 = app.Ei0EditField.Value;
            w = app.EditField.Value;
            theta_i = app.iEditField.Value;
            epsilon1 = app.r1EditField.Value;
            mu1 = app.r1EditField_2.Value;
            epsilon2 = app.r2EditField.Value;
            mu2 = app.r2EditField_2.Value;
            
            eta1 = 120*pi*sqrt(mu1/epsilon1);
            %eta1
            eta2 = 120*pi*sqrt(mu2/epsilon2);
            %eta2
            e0 = 8.8541878128*10^-12;
            m0 = 4*pi*10^-7;
            beta1 = w*sqrt(mu1*epsilon1*e0*m0);
            %beta1
            beta2 = w*sqrt(mu2*epsilon2*e0*m0);
            %beta2
            theta_r = theta_i;
            %theta_r
            theta_t = asind(sind(theta_i)*beta1/beta2);
            %theta_t
            gamma1 = ((eta2*cosd(theta_t)) - (eta1*cosd(theta_i)))/((eta2*cosd(theta_t)) + (eta1*cosd(theta_i)));
            %gamma1
            x = linspace(-0.5,0,101);
            z = linspace(-0.5,0,101);
            y = linspace(-0.5,0,101);
            %defining E field Variable
            Eix = zeros(101,101);
            Eiy = zeros(101,101);
            Eiz = zeros(101,101);
            Erx = zeros(101,101);
            Ery = zeros(101,101);
            Erz = zeros(101,101);
            %setting values
            for n = 1:1:101
                for m = 1:1:101
                    Eix(n, m) = sind(theta_i).*Ei0.*exp(-1i*(beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                    Eiz(n, m) = -cosd(theta_i).*Ei0.*exp(-1i*(beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                    Erx(n, m) = sind(theta_i).*gamma1.*Ei0.*exp(-1i*(-beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                    Erz(n, m) = cosd(theta_i).*gamma1.*Ei0.*exp(-1i*(-beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                end
            end
            subplot(2,2,1);
            surf(x, z, real(Eix + Eiz));
            title("Parallel Polarization electric field")
            subtitle('Ei ');
            xlabel("X-axis")
            ylabel("Z-axis")
            zlabel("Y-axis")
            subplot(2,2,2);
            surf(x, z, real(Erx + Erz));
            
            subtitle('Er');
            xlabel("X-axis")
            ylabel("Z-axis")
            zlabel("Y-axis")
            subplot(2,2,[3,4]);
            surf(x, z, real(Eix + Erx + Eiz + Erz));
            
            subtitle('Total E1(Ei + Er)');
            xlabel("X-axis")
            ylabel("Z-axis")
            zlabel("Y-axis")
        end

        % Button pushed function: PlotE2EtButton_2
        function PlotE2EtButton_2Pushed(app, event)
            Ei0 = app.Ei0EditField.Value;
            w = app.EditField.Value;
            theta_i = app.iEditField.Value;
            epsilon1 = app.r1EditField.Value;
            mu1 = app.r1EditField_2.Value;
            epsilon2 = app.r2EditField.Value;
            mu2 = app.r2EditField_2.Value;
            
            eta1 = 120*pi*sqrt(mu1/epsilon1);
            %eta1
            eta2 = 120*pi*sqrt(mu2/epsilon2);
            %eta2
            ep = 8.8541878128*10^-12;
            m = 4*pi*10^-7;
            beta1 = w*sqrt(mu1*epsilon1*ep*m);
            %beta1
            beta2 = w*sqrt(mu2*epsilon2*ep*m);
            %beta2
            
            theta_t = asind(sind(theta_i)*beta1/beta2);
            %theta_t
            
            tau1 = (2*eta2*cosd(theta_i))/((eta2*cosd(theta_t)) + (eta1*cosd(theta_i)));
            %tau1
            
            %Plot for E2
            x = linspace(0,0.5,101);
            z = linspace(0,0.5,101);
            y = linspace(0,0.5,101);
            %defining E field Variable
            Etx = zeros(101,101);
            Ety = zeros(101,101);
            Etz = zeros(101,101);
            
            %setting values
            for n = 1:1:101
                for m = 1:1:101
                    Etx(n, m) = sind(theta_t).*tau1.*Ei0.*exp(-1i*(beta2.*cosd(theta_t).*x(n) + beta2.*sind(theta_t).*z(m)));
                    Etz(n, m) = -cosd(theta_t).*tau1.*Ei0.*exp(-1i*(beta2.*cosd(theta_t).*x(n) + beta2.*sind(theta_t).*z(m)));
                    
                end
            end
            surf(x,z,real(Etx+Etz));
            title('Parallel Polarization Electric field');
            subtitle('E2(Et)');
            xlabel('X-axis');
            ylabel('Z-axis');
            zlabel('Y-axis');
        end

        % Button pushed function: PlotH1HiHrButton_2
        function PlotH1HiHrButton_2Pushed(app, event)
            Ei0 = app.Ei0EditField.Value;
            w = app.EditField.Value;
            theta_i = app.iEditField.Value;
            epsilon1 = app.r1EditField.Value;
            mu1 = app.r1EditField_2.Value;
            epsilon2 = app.r2EditField.Value;
            mu2 = app.r2EditField_2.Value;
            
            eta1 = 120*pi*sqrt(mu1/epsilon1);
            %eta1
            eta2 = 120*pi*sqrt(mu2/epsilon2);
            %eta2
            e0 = 8.8541878128*10^-12;
            m0 = 4*pi*10^-7;
            beta1 = w*sqrt(mu1*epsilon1*e0*m0);
            %beta1
            beta2 = w*sqrt(mu2*epsilon2*e0*m0);
            %beta2
            theta_r = theta_i;
            %theta_r
            theta_t = asind(sind(theta_i)*beta1/beta2);
            %theta_t
            gamma1 = ((eta2*cosd(theta_t)) - (eta1*cosd(theta_i)))/((eta2*cosd(theta_t)) + (eta1*cosd(theta_i)));
            %gamma1
            x = linspace(-0.5,0,101);
            z = linspace(-0.5,0,101);
            y = linspace(-0.5,0,101);
            Hix = zeros(101, 101);
            Hiy = zeros(101, 101);
            Hiz = zeros(101, 101);
            
            Hrx = zeros(101, 101);
            Hry = zeros(101, 101);
            Hrz = zeros(101, 101);
            for n = 1:1:101
                for m = 1:1:101
                    Hiy(n, m) = (Ei0/eta1).*exp(-1i*(beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                    Hry(n, m) = (gamma1.*Ei0/eta1).*exp(-1i*(-beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                end
            end
            subplot(2,2,1);
            surf(x, z, real(Hiy));
            title('Parallel Polarization Magnetic Field');
            subtitle('Hi');
            xlabel('X-axis');
            ylabel('Z-axis');
            zlabel('Y-axis');
            subplot(2,2,2);
            surf(x, z, real(Hry));
            
            subtitle('Hr');
            xlabel('X-axis');
            ylabel('Z-axis');
            zlabel('Y-axis');
            subplot(2,2,[3,4]);
            surf(x, z, real(Hiy + Hry ));
            
            subtitle('Total H1(Hi + Hr)');
            xlabel("X-axis")
            ylabel("Z-axis")
            zlabel("Y-axis")
        end

        % Button pushed function: PlotH2HtButton_2
        function PlotH2HtButton_2Pushed(app, event)
            Ei0 = app.Ei0EditField.Value;
            w = app.EditField.Value;
            theta_i = app.iEditField.Value;
            epsilon1 = app.r1EditField.Value;
            mu1 = app.r1EditField_2.Value;
            epsilon2 = app.r2EditField.Value;
            mu2 = app.r2EditField_2.Value;
            
            eta1 = 120*pi*sqrt(mu1/epsilon1);
            %eta1
            eta2 = 120*pi*sqrt(mu2/epsilon2);
            %eta2
            ep = 8.8541878128*10^-12;
            m = 4*pi*10^-7;
            beta1 = w*sqrt(mu1*epsilon1*ep*m);
            %beta1
            beta2 = w*sqrt(mu2*epsilon2*ep*m);
            %beta2
            
            theta_t = asind(sind(theta_i)*beta1/beta2);
            %theta_t
            
            tau1 = (2*eta2*cosd(theta_i))/((eta2*cosd(theta_t)) + (eta1*cosd(theta_i)));
            %tau1
            
            %Plot for E2
            x = linspace(0,0.5,101);
            z = linspace(0,0.5,101);
            y = linspace(0,0.5,101);
            %defining E field Variable
            Htx = zeros(101,101);
            Hty = zeros(101,101);
            Htz = zeros(101,101);
            
            %setting values
            for n = 1:1:101
                for m = 1:1:101
                    Hty(n, m) = (tau1.*Ei0/eta2).*exp(-1i*(beta2.*cosd(theta_t).*x(n) + beta2.*sind(theta_t).*z(m)));
                end
            end
            surf(x, z, real(Hty));
            title('Parallel Polarization Magnetic Field');
            subtitle('Total H2(Ht)');
            xlabel('X-axis');
            ylabel('Z-axis');
            zlabel('Y-axis');
        end

        % Button pushed function: PlotE1E2Button
        function PlotE1E2ButtonPushed(app, event)
            Ei0 = app.Ei0EditField.Value;
            w = app.EditField.Value;
            theta_i = app.iEditField.Value;
            epsilon1 = app.r1EditField.Value;
            mu1 = app.r1EditField_2.Value;
            epsilon2 = app.r2EditField.Value;
            mu2 = app.r2EditField_2.Value;
            
            eta1 = 120*pi*sqrt(mu1/epsilon1);
            %eta1
            eta2 = 120*pi*sqrt(mu2/epsilon2);
            %eta2
            e0 = 8.8541878128*10^-12;
            m0 = 4*pi*10^-7;
            beta1 = w*sqrt(mu1*epsilon1*e0*m0);
            %beta1
            beta2 = w*sqrt(mu2*epsilon2*e0*m0);
            %beta2
            theta_r = theta_i;
            %theta_r
            theta_t = asind(sind(theta_i)*beta1/beta2);
            %theta_t
            gamma1 = ((eta2*cosd(theta_i)) - (eta1*cosd(theta_t)))/((eta2*cosd(theta_i)) + (eta1*cosd(theta_t)));
            %gamma1
            tau1 = (2*eta2*cosd(theta_i))/((eta2*cosd(theta_i)) + (eta1*cosd(theta_t)));
            %tau1
             x = linspace(-0.5,0,101);
            z = linspace(-0.5,0,101);
            y = linspace(-0.5,0,101);
            %defining E field Variable
            Eix = zeros(101,101);
            Eiy = zeros(101,101);
            Eiz = zeros(101,101);
            Erx = zeros(101,101);
            Ery = zeros(101,101);
            Erz = zeros(101,101);
            %setting values
            for n = 1:1:101
                for m = 1:1:101
                    Eiy(n, m) = Ei0.*exp(-1i*(beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                    Ery(n, m) = Ei0.*gamma1.*exp(-1i*(-beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                end
            end
            x1 = linspace(0,0.5,101);
            z1 = linspace(0,0.5,101);
            y1 = linspace(0,0.5,101);
            %defining E field Variable
            Etx = zeros(101,101);
            Ety = zeros(101,101);
            Etz = zeros(101,101);
            
            %setting values
            for n = 1:1:101
                for m = 1:1:101
                    Ety(n, m) = tau1.*Ei0.*exp(-1i*(beta2.*cosd(theta_t).*x1(n) + beta2.*sind(theta_t).*z1(m)));
                    
                end
            end
            subplot(2,2,[1,2]);
            surf(x,z,real(Eiy+Ery));
            title('Perpendicular Polarization Electric Field');
            subtitle('E1');
            xlabel('X-axis');
            ylabel('Z-axis');
            zlabel('Y-axis');
            subplot(2,2,[3,4]);
            surf(x1,z1,real(Ety));
            
            subtitle('E2');
            xlabel('X-axis');
            ylabel('Z-axis');
            zlabel('Y-axis');
        end

        % Button pushed function: PlotH1H2Button
        function PlotH1H2ButtonPushed(app, event)
            Ei0 = app.Ei0EditField.Value;
            w = app.EditField.Value;
            theta_i = app.iEditField.Value;
            epsilon1 = app.r1EditField.Value;
            mu1 = app.r1EditField_2.Value;
            epsilon2 = app.r2EditField.Value;
            mu2 = app.r2EditField_2.Value;
            
            eta1 = 120*pi*sqrt(mu1/epsilon1);
            %eta1
            eta2 = 120*pi*sqrt(mu2/epsilon2);
            %eta2
            e0 = 8.8541878128*10^-12;
            m0 = 4*pi*10^-7;
            beta1 = w*sqrt(mu1*epsilon1*e0*m0);
            %beta1
            beta2 = w*sqrt(mu2*epsilon2*e0*m0);
            %beta2
            theta_r = theta_i;
            %theta_r
            theta_t = asind(sind(theta_i)*beta1/beta2);
            %theta_t
            gamma1 = ((eta2*cosd(theta_i)) - (eta1*cosd(theta_t)))/((eta2*cosd(theta_i)) + (eta1*cosd(theta_t)));
            %gamma1
            tau1 = (2*eta2*cosd(theta_i))/((eta2*cosd(theta_i)) + (eta1*cosd(theta_t)));
            %tau1
            x = linspace(-0.5,0,101);
            z = linspace(-0.5,0,101);
            y = linspace(-0.5,0,101);
            Hix = zeros(101, 101);
            Hiy = zeros(101, 101);
            Hiz = zeros(101, 101);
            
            Hrx = zeros(101, 101);
            Hry = zeros(101, 101);
            Hrz = zeros(101, 101);
            for n = 1:1:101
                for m = 1:1:101
                    Hix(n, m) = -sind(theta_i).*(Ei0/eta1).*exp(-1i*(beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                    Hiz(n, m) = cosd(theta_i).*(Ei0/eta1).*exp(-1i*(beta1.*cosd(theta_i).*x(n) +  beta1.*sind(theta_i).*z(m)));
                    Hrx(n, m) = -sind(theta_i).*(gamma1.*Ei0/eta1).*exp(-1i*(-beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                    Hrz(n, m) = -cosd(theta_i).*(gamma1.*Ei0/eta1).*exp(-1i*(-beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                end
            end
            %Plot for E2
            x1 = linspace(0,0.5,101);
            z1 = linspace(0,0.5,101);
            y1 = linspace(0,0.5,101);
            %defining E field Variable
            Htx = zeros(101,101);
            Hty = zeros(101,101);
            Htz = zeros(101,101);
            
            %setting values
            for n = 1:1:101
                for m = 1:1:101
                    Htx(n, m) = -sind(theta_t).*(tau1.*Ei0/eta2).*exp(-1i*(beta2.*cosd(theta_t).*x1(n) + beta2.*sind(theta_t).*z1(m)));
                    Htz(n, m) = cosd(theta_t).*(tau1.*Ei0/eta2).*exp(-1i*(beta2.*cosd(theta_t).*x1(n) +  beta2.*sind(theta_t).*z1(m)));
                end
            end
            
            subplot(2,2,[1,2]);
            surf(x,z,real(Hix + Hiz + Hrx + Hrz));
            title('Perpendicular Polarization Magnetic Field');
            subtitle('H1');
            xlabel('X-axis');
            ylabel('Z-axis');
            zlabel('Y-axis');
            subplot(2,2,[3,4]);
            surf(x1,z1,real(Htx + Htz));
            subtitle('H2');
            xlabel('X-axis');
            ylabel('Z-axis');
            zlabel('Y-axis');
            
            
        end

        % Button pushed function: PlotE1E2Button_2
        function PlotE1E2Button_2Pushed(app, event)
            Ei0 = app.Ei0EditField.Value;
            w = app.EditField.Value;
            theta_i = app.iEditField.Value;
            epsilon1 = app.r1EditField.Value;
            mu1 = app.r1EditField_2.Value;
            epsilon2 = app.r2EditField.Value;
            mu2 = app.r2EditField_2.Value;
            
            eta1 = 120*pi*sqrt(mu1/epsilon1);
            %eta1
            eta2 = 120*pi*sqrt(mu2/epsilon2);
            %eta2
            e0 = 8.8541878128*10^-12;
            m0 = 4*pi*10^-7;
            beta1 = w*sqrt(mu1*epsilon1*e0*m0);
            %beta1
            beta2 = w*sqrt(mu2*epsilon2*e0*m0);
            %beta2
            theta_r = theta_i;
            %theta_r
            theta_t = asind(sind(theta_i)*beta1/beta2);
            %theta_t
            gamma1 = ((eta2*cosd(theta_t)) - (eta1*cosd(theta_i)))/((eta2*cosd(theta_t)) + (eta1*cosd(theta_i)));
            %gamma1
            tau1 = (2*eta2*cosd(theta_i))/((eta2*cosd(theta_t)) + (eta1*cosd(theta_i)));
            %tau1
            x = linspace(-0.5,0,101);
            z = linspace(-0.5,0,101);
            y = linspace(-0.5,0,101);
            %defining E field Variable
            Eix = zeros(101,101);
            Eiy = zeros(101,101);
            Eiz = zeros(101,101);
            Erx = zeros(101,101);
            Ery = zeros(101,101);
            Erz = zeros(101,101);
            %setting values
            for n = 1:1:101
                for m = 1:1:101
                    Eix(n, m) = sind(theta_i).*Ei0.*exp(-1i*(beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                    Eiz(n, m) = -cosd(theta_i).*Ei0.*exp(-1i*(beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                    Erx(n, m) = sind(theta_i).*gamma1.*Ei0.*exp(-1i*(-beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                    Erz(n, m) = cosd(theta_i).*gamma1.*Ei0.*exp(-1i*(-beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                end
            end
            
            
            %Plot for E2
            x1 = linspace(0,0.5,101);
            z1 = linspace(0,0.5,101);
            y1 = linspace(0,0.5,101);
            %defining E field Variable
            Etx = zeros(101,101);
            Ety = zeros(101,101);
            Etz = zeros(101,101);
            
            %setting values
            for n = 1:1:101
                for m = 1:1:101
                    Etx(n, m) = sind(theta_t).*tau1.*Ei0.*exp(-1i*(beta2.*cosd(theta_t).*x1(n) + beta2.*sind(theta_t).*z1(m)));
                    Etz(n, m) = -cosd(theta_t).*tau1.*Ei0.*exp(-1i*(beta2.*cosd(theta_t).*x1(n) + beta2.*sind(theta_t).*z1(m)));
                    
                end
            end
            subplot(2,2,[1,2]);
            surf(x, z, real(Eix + Erx + Eiz + Erz));
            title('Parallel Polarization Electric Field')
            subtitle('E1');
            xlabel("X-axis")
            ylabel("Z-axis")
            zlabel("Y-axis")
            subplot(2,2,[3,4]);
            surf(x, z, real(Etx + Etz));
            subtitle('E2');
            xlabel("X-axis")
            ylabel("Z-axis")
            zlabel("Y-axis")
        end

        % Button pushed function: PlotH1H2Button_2
        function PlotH1H2Button_2Pushed(app, event)
            Ei0 = app.Ei0EditField.Value;
            w = app.EditField.Value;
            theta_i = app.iEditField.Value;
            epsilon1 = app.r1EditField.Value;
            mu1 = app.r1EditField_2.Value;
            epsilon2 = app.r2EditField.Value;
            mu2 = app.r2EditField_2.Value;
            
            eta1 = 120*pi*sqrt(mu1/epsilon1);
            %eta1
            eta2 = 120*pi*sqrt(mu2/epsilon2);
            %eta2
            e0 = 8.8541878128*10^-12;
            m0 = 4*pi*10^-7;
            beta1 = w*sqrt(mu1*epsilon1*e0*m0);
            %beta1
            beta2 = w*sqrt(mu2*epsilon2*e0*m0);
            %beta2
            theta_r = theta_i;
            %theta_r
            theta_t = asind(sind(theta_i)*beta1/beta2);
            %theta_t
            gamma1 = ((eta2*cosd(theta_t)) - (eta1*cosd(theta_i)))/((eta2*cosd(theta_t)) + (eta1*cosd(theta_i)));
            %gamma1
            tau1 = (2*eta2*cosd(theta_i))/((eta2*cosd(theta_t)) + (eta1*cosd(theta_i)));
            %tau1
            x = linspace(-0.5,0,101);
            z = linspace(-0.5,0,101);
            y = linspace(-0.5,0,101);
            Hix = zeros(101, 101);
            Hiy = zeros(101, 101);
            Hiz = zeros(101, 101);
            
            Hrx = zeros(101, 101);
            Hry = zeros(101, 101);
            Hrz = zeros(101, 101);
            for n = 1:1:101
                for m = 1:1:101
                    Hiy(n, m) = (Ei0/eta1).*exp(-1i*(beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                    Hry(n, m) = (gamma1.*Ei0/eta1).*exp(-1i*(-beta1.*cosd(theta_i).*x(n) + beta1.*sind(theta_i).*z(m)));
                end
            end
            
            %Plot for E2
            x1 = linspace(0,0.5,101);
            z1 = linspace(0,0.5,101);
            y1 = linspace(0,0.5,101);
            %defining E field Variable
            Htx = zeros(101,101);
            Hty = zeros(101,101);
            Htz = zeros(101,101);
            
            %setting values
            for n = 1:1:101
                for m = 1:1:101
                    Hty(n, m) = (tau1.*Ei0/eta2).*exp(-1i*(beta2.*cosd(theta_t).*x1(n) + beta2.*sind(theta_t).*z1(m)));
                end
            end
            subplot(2,2,[1,2]);
            surf(x, z, real(Hiy + Hry));
            title('Parallel Polarization Magnetic Field')
            subtitle('H1');
            xlabel("X-axis")
            ylabel("Z-axis")
            zlabel("Y-axis")
            subplot(2,2,[3,4]);
            surf(x, z, real(Hty));
            subtitle('H2');
            xlabel("X-axis")
            ylabel("Z-axis")
            zlabel("Y-axis")
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.8 0.8 0.8];
            app.UIFigure.Position = [100 100 646 536];
            app.UIFigure.Name = 'MATLAB App';

            % Create ObliqueIncidenceLabel
            app.ObliqueIncidenceLabel = uilabel(app.UIFigure);
            app.ObliqueIncidenceLabel.HorizontalAlignment = 'center';
            app.ObliqueIncidenceLabel.WordWrap = 'on';
            app.ObliqueIncidenceLabel.FontSize = 16;
            app.ObliqueIncidenceLabel.FontWeight = 'bold';
            app.ObliqueIncidenceLabel.Position = [255 494 150 37];
            app.ObliqueIncidenceLabel.Text = 'Oblique Incidence';

            % Create WaveParametersLabel
            app.WaveParametersLabel = uilabel(app.UIFigure);
            app.WaveParametersLabel.FontWeight = 'bold';
            app.WaveParametersLabel.Position = [68 473 105 22];
            app.WaveParametersLabel.Text = 'Wave Parameters';

            % Create Ei0EditFieldLabel
            app.Ei0EditFieldLabel = uilabel(app.UIFigure);
            app.Ei0EditFieldLabel.HorizontalAlignment = 'right';
            app.Ei0EditFieldLabel.Position = [50 432 25 22];
            app.Ei0EditFieldLabel.Text = 'Ei0';

            % Create Ei0EditField
            app.Ei0EditField = uieditfield(app.UIFigure, 'numeric');
            app.Ei0EditField.Position = [90 432 100 22];

            % Create EditFieldLabel
            app.EditFieldLabel = uilabel(app.UIFigure);
            app.EditFieldLabel.HorizontalAlignment = 'right';
            app.EditFieldLabel.Position = [51 395 25 22];
            app.EditFieldLabel.Text = 'ω';

            % Create EditField
            app.EditField = uieditfield(app.UIFigure, 'numeric');
            app.EditField.Position = [91 395 100 22];

            % Create iEditFieldLabel
            app.iEditFieldLabel = uilabel(app.UIFigure);
            app.iEditFieldLabel.HorizontalAlignment = 'right';
            app.iEditFieldLabel.Position = [51 356 25 22];
            app.iEditFieldLabel.Text = 'θi';

            % Create iEditField
            app.iEditField = uieditfield(app.UIFigure, 'numeric');
            app.iEditField.Position = [91 356 100 22];

            % Create r1EditFieldLabel
            app.r1EditFieldLabel = uilabel(app.UIFigure);
            app.r1EditFieldLabel.HorizontalAlignment = 'right';
            app.r1EditFieldLabel.Position = [255 432 25 22];
            app.r1EditFieldLabel.Text = 'εr1';

            % Create r1EditField
            app.r1EditField = uieditfield(app.UIFigure, 'numeric');
            app.r1EditField.Position = [295 432 100 22];

            % Create r1EditField_2Label
            app.r1EditField_2Label = uilabel(app.UIFigure);
            app.r1EditField_2Label.HorizontalAlignment = 'right';
            app.r1EditField_2Label.Position = [483 432 25 22];
            app.r1EditField_2Label.Text = 'μr1';

            % Create r1EditField_2
            app.r1EditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.r1EditField_2.Position = [523 432 100 22];

            % Create r2EditFieldLabel
            app.r2EditFieldLabel = uilabel(app.UIFigure);
            app.r2EditFieldLabel.HorizontalAlignment = 'right';
            app.r2EditFieldLabel.Position = [255 395 25 22];
            app.r2EditFieldLabel.Text = 'εr2';

            % Create r2EditField
            app.r2EditField = uieditfield(app.UIFigure, 'numeric');
            app.r2EditField.Position = [295 395 100 22];

            % Create r2EditField_2Label
            app.r2EditField_2Label = uilabel(app.UIFigure);
            app.r2EditField_2Label.HorizontalAlignment = 'right';
            app.r2EditField_2Label.Position = [484 395 25 22];
            app.r2EditField_2Label.Text = 'μr2';

            % Create r2EditField_2
            app.r2EditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.r2EditField_2.Position = [524 395 100 22];

            % Create PerpendicularPolarizationLabel
            app.PerpendicularPolarizationLabel = uilabel(app.UIFigure);
            app.PerpendicularPolarizationLabel.HorizontalAlignment = 'center';
            app.PerpendicularPolarizationLabel.FontSize = 13;
            app.PerpendicularPolarizationLabel.FontWeight = 'bold';
            app.PerpendicularPolarizationLabel.Position = [85 317 173 40];
            app.PerpendicularPolarizationLabel.Text = 'Perpendicular  Polarization';

            % Create PlotE1EiErButton
            app.PlotE1EiErButton = uibutton(app.UIFigure, 'push');
            app.PlotE1EiErButton.ButtonPushedFcn = createCallbackFcn(app, @PlotE1EiErButtonPushed, true);
            app.PlotE1EiErButton.BackgroundColor = [0.9412 0.9412 0.9412];
            app.PlotE1EiErButton.Position = [111 295 100 22];
            app.PlotE1EiErButton.Text = 'Plot E1(Ei + Er)';

            % Create PlotE2EtButton
            app.PlotE2EtButton = uibutton(app.UIFigure, 'push');
            app.PlotE2EtButton.ButtonPushedFcn = createCallbackFcn(app, @PlotE2EtButtonPushed, true);
            app.PlotE2EtButton.Position = [111 265 100 22];
            app.PlotE2EtButton.Text = 'Plot E2(Et)';

            % Create PlotH1HiHrButton
            app.PlotH1HiHrButton = uibutton(app.UIFigure, 'push');
            app.PlotH1HiHrButton.ButtonPushedFcn = createCallbackFcn(app, @PlotH1HiHrButtonPushed, true);
            app.PlotH1HiHrButton.Position = [111 235 101 22];
            app.PlotH1HiHrButton.Text = 'Plot H1(Hi + Hr)';

            % Create PlotH2HtButton
            app.PlotH2HtButton = uibutton(app.UIFigure, 'push');
            app.PlotH2HtButton.ButtonPushedFcn = createCallbackFcn(app, @PlotH2HtButtonPushed, true);
            app.PlotH2HtButton.Position = [111 205 100 22];
            app.PlotH2HtButton.Text = 'Plot H2(Ht)';

            % Create BrewsterAngleBEditFieldLabel
            app.BrewsterAngleBEditFieldLabel = uilabel(app.UIFigure);
            app.BrewsterAngleBEditFieldLabel.HorizontalAlignment = 'right';
            app.BrewsterAngleBEditFieldLabel.Position = [28 70 110 22];
            app.BrewsterAngleBEditFieldLabel.Text = 'Brewster Angle(θB)';

            % Create BrewsterAngleBEditField
            app.BrewsterAngleBEditField = uieditfield(app.UIFigure, 'numeric');
            app.BrewsterAngleBEditField.Position = [153 70 100 22];

            % Create ParallelPolarizationLabel
            app.ParallelPolarizationLabel = uilabel(app.UIFigure);
            app.ParallelPolarizationLabel.HorizontalAlignment = 'center';
            app.ParallelPolarizationLabel.FontSize = 13;
            app.ParallelPolarizationLabel.FontWeight = 'bold';
            app.ParallelPolarizationLabel.Position = [411 317 160 40];
            app.ParallelPolarizationLabel.Text = 'Parallel Polarization';

            % Create PlotE1EiErButton_2
            app.PlotE1EiErButton_2 = uibutton(app.UIFigure, 'push');
            app.PlotE1EiErButton_2.ButtonPushedFcn = createCallbackFcn(app, @PlotE1EiErButton_2Pushed, true);
            app.PlotE1EiErButton_2.Position = [441 295 100 22];
            app.PlotE1EiErButton_2.Text = 'Plot E1(Ei + Er)';

            % Create PlotE2EtButton_2
            app.PlotE2EtButton_2 = uibutton(app.UIFigure, 'push');
            app.PlotE2EtButton_2.ButtonPushedFcn = createCallbackFcn(app, @PlotE2EtButton_2Pushed, true);
            app.PlotE2EtButton_2.Position = [441 265 100 22];
            app.PlotE2EtButton_2.Text = 'Plot E2(Et)';

            % Create PlotH1HiHrButton_2
            app.PlotH1HiHrButton_2 = uibutton(app.UIFigure, 'push');
            app.PlotH1HiHrButton_2.ButtonPushedFcn = createCallbackFcn(app, @PlotH1HiHrButton_2Pushed, true);
            app.PlotH1HiHrButton_2.Position = [441 235 101 22];
            app.PlotH1HiHrButton_2.Text = 'Plot H1(Hi + Hr)';

            % Create PlotH2HtButton_2
            app.PlotH2HtButton_2 = uibutton(app.UIFigure, 'push');
            app.PlotH2HtButton_2.ButtonPushedFcn = createCallbackFcn(app, @PlotH2HtButton_2Pushed, true);
            app.PlotH2HtButton_2.Position = [441 205 100 22];
            app.PlotH2HtButton_2.Text = 'Plot H2(Ht)';

            % Create CalculateBrewsterAngleButton
            app.CalculateBrewsterAngleButton = uibutton(app.UIFigure, 'push');
            app.CalculateBrewsterAngleButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateBrewsterAngleButtonPushed, true);
            app.CalculateBrewsterAngleButton.Position = [90 102 151 22];
            app.CalculateBrewsterAngleButton.Text = 'Calculate Brewster Angle';

            % Create CalculateBrewsterAngleButton_2
            app.CalculateBrewsterAngleButton_2 = uibutton(app.UIFigure, 'push');
            app.CalculateBrewsterAngleButton_2.ButtonPushedFcn = createCallbackFcn(app, @CalculateBrewsterAngleButton_2Pushed, true);
            app.CalculateBrewsterAngleButton_2.Position = [416 102 151 22];
            app.CalculateBrewsterAngleButton_2.Text = 'Calculate Brewster Angle';

            % Create BrewsterAngleBEditField_2Label
            app.BrewsterAngleBEditField_2Label = uilabel(app.UIFigure);
            app.BrewsterAngleBEditField_2Label.HorizontalAlignment = 'right';
            app.BrewsterAngleBEditField_2Label.Position = [346 70 110 22];
            app.BrewsterAngleBEditField_2Label.Text = 'Brewster Angle(θB)';

            % Create BrewsterAngleBEditField_2
            app.BrewsterAngleBEditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.BrewsterAngleBEditField_2.Position = [471 70 100 22];

            % Create degreeLabel
            app.degreeLabel = uilabel(app.UIFigure);
            app.degreeLabel.Position = [262 70 43 22];
            app.degreeLabel.Text = 'degree';

            % Create degreeLabel_2
            app.degreeLabel_2 = uilabel(app.UIFigure);
            app.degreeLabel_2.Position = [580 70 43 22];
            app.degreeLabel_2.Text = 'degree';

            % Create degreeLabel_3
            app.degreeLabel_3 = uilabel(app.UIFigure);
            app.degreeLabel_3.Position = [208 356 43 22];
            app.degreeLabel_3.Text = 'degree';

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.HorizontalAlignment = 'center';
            app.Label.FontWeight = 'bold';
            app.Label.FontColor = [1 0 0];
            app.Label.Position = [68 20 527 35];
            app.Label.Text = '-';

            % Create MediumParametersLabel
            app.MediumParametersLabel = uilabel(app.UIFigure);
            app.MediumParametersLabel.FontWeight = 'bold';
            app.MediumParametersLabel.Position = [393 473 120 22];
            app.MediumParametersLabel.Text = 'Medium Parameters';

            % Create PlotE1E2Button
            app.PlotE1E2Button = uibutton(app.UIFigure, 'push');
            app.PlotE1E2Button.ButtonPushedFcn = createCallbackFcn(app, @PlotE1E2ButtonPushed, true);
            app.PlotE1E2Button.Position = [112 170 100 22];
            app.PlotE1E2Button.Text = 'Plot E1, E2';

            % Create PlotE1E2Button_2
            app.PlotE1E2Button_2 = uibutton(app.UIFigure, 'push');
            app.PlotE1E2Button_2.ButtonPushedFcn = createCallbackFcn(app, @PlotE1E2Button_2Pushed, true);
            app.PlotE1E2Button_2.Position = [442 170 100 22];
            app.PlotE1E2Button_2.Text = 'Plot E1, E2';

            % Create PlotH1H2Button
            app.PlotH1H2Button = uibutton(app.UIFigure, 'push');
            app.PlotH1H2Button.ButtonPushedFcn = createCallbackFcn(app, @PlotH1H2ButtonPushed, true);
            app.PlotH1H2Button.Position = [112 134 100 22];
            app.PlotH1H2Button.Text = 'Plot H1, H2';

            % Create PlotH1H2Button_2
            app.PlotH1H2Button_2 = uibutton(app.UIFigure, 'push');
            app.PlotH1H2Button_2.ButtonPushedFcn = createCallbackFcn(app, @PlotH1H2Button_2Pushed, true);
            app.PlotH1H2Button_2.Position = [446 134 100 22];
            app.PlotH1H2Button_2.Text = 'Plot H1, H2';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Project_12_final_draft_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end