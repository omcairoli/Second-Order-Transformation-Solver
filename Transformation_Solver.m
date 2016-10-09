%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%              Second-Order Transformation Solver               %%%%
%%%%          Just run the program...Yes, it's that easy           %%%%
%%%%                         TOP SECRET                            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;clc;
close all;
format compact;
format short;

x = cell2mat(inputdlg('Enter coefficients of Numerator:'));
num = sscanf(x,'%f %f %f');
num = num';

y = cell2mat(inputdlg('Enter coefficients of Denominator:'));
den = sscanf(y,'%f %f %f');
den = den';

z = str2double(cell2mat(inputdlg('Enter sampling frequency (fsampling):')));

% ----------------------------------------------------------------------- %
fs = z;                                       %% Sampling Frequency
d = num;                                      %% Numerator coefficients
c = den;                                      %% Denominator coefficients
% ------------------------------------------------------------------------%


Approx_List = {'Impulse Invariant Transformation','Bilinear Transformation','Euler Approximation'};
Approx_Select = listdlg('ListString',Approx_List','SelectionMode','single','ListSize',[300 75],'Name','Select Transformation:');

if (Approx_Select == 1)
    
    HPF = 0;
    BSF = 0;
    
    if (d(1)>0)
        d_1 = [0 0 1];
        c_1 = c;
        N1_1 = length(c_1);
        N2_1 = length(d_1);
        d_1 = [zeros(1,N1_1-N2_1),d_1];
        [R_1,AP_1,C_1] = residue(d_1,c_1);
        [d1_1,c1_1] = residue(R_1,AP_1,C_1); % checks
        AZ_1 = roots(d1_1);
        T_1 = 1/fs;
        DP_1 = exp(AP_1*T_1);
        [b_1,a_1] = residuez(R_1,DP_1,C_1);
        [R1_1,P1_1,C1_1] = residuez(b_1,a_1);
        b_1 = [b_1,zeros(1,length(a_1)-length(b_1))];
        DZ_1 = roots(b_1);
        
        d_2 = [0 0 1];
        c_2 = d;
        N1_2 = length(c_2);
        N2_2 = length(d_2);
        d_2 = [zeros(1,N1_2-N2_2),d_2];
        [R_2,AP_2,C_2] = residue(d_2,c_2);
        [d1_2,c1_2] = residue(R_2,AP_2,C_2); % checks
        AZ_2 = roots(d1_2);
        T_2 = 1/fs;
        DP_2 = exp(AP_2*T_2);
        [b_2,a_2] = residuez(R_2,DP_2,C_2);
        [R1_2,P1_2,C1_2] = residuez(b_2,a_2);
        b_2 = [b_2,zeros(1,length(a_2)-length(b_2))];
        DZ_2 = roots(b_2);
        
        if (d(3)==0)
            HPF = 1;
        end
        if (d(3)>0)
            BSF = 1;
        end
        
        AP = AP_1;
        AZ = AP_2;
        N1 = N1_1;
        T = T_1;
        
    end
    
    if (BSF==0 && HPF==0)
        N1 = length(c);
        N2 = length(d);
        d = [zeros(1,N1-N2),d];
        [R,AP,C] = residue(d,c);
        [d1,c1] = residue(R,AP,C); % checks
        AZ = roots(d1);
        T = 1/fs;
        DP = exp(AP*T);
        [b,a] = residuez(R,DP,C);
        [R1,P1,C1] = residuez(b,a);
        b = [b,zeros(1,length(a)-length(b))];
        DZ = roots(b);
    end
    
    %% Analog Pole-Zero Diagram
    
    figure(1);
    plot(real(AP),imag(AP),'bx','MarkerSize',12);hold on;
    plot(real(AZ),imag(AZ),'bo','MarkerSize',12);grid on;
    title('Analog Pole-Zero Diagram');
    xlabel('Real Part');
    ylabel('Imaginary Part');
    axis square;
    
    %% Analog Magnitude and Phase Response
    
    wmax = 2*pi*fs/2;
    w = 0:wmax/1000:wmax;
    E = (N1-1:-1:0)'*ones(1,length(w));
    W = ones(1,N1)'*1i*w;
    W = W.^E;
    Ha = d*W./(c*W);
    w1 = logspace(-2,2,length(Ha));
    E1 = (N1-1:-1:0)'*ones(1,length(w1));
    W1 = ones(1,N1)'*1i*w1;
    W1 = W1.^E1;
    Hb = d*W1./(c*W1);
    figure(2);
    plot(w,abs(Ha));grid on;
    xlabel('\omega');ylabel('|H(\omega)|');
    title('Analog Magnitude Response (Linear)');
    figure(3);
    semilogx(w1,20*log10(abs(Hb)));grid on;
    xlabel('\omega');ylabel('20*log10|H(\omega)|');
    title('Analog Magnitude Response (dB)');
    figure(4);
    plot(w,180/pi*angle(Ha));grid on;
    xlabel('\omega');ylabel('Angle(H(\omega) (deg)');
    title('Analog Phase Response');
    
    
    %% Digital Pole-Zero Diagram
    
    x = -1:0.01:1;
    ty = sqrt(1-x.^2);
    figure(5);
    plot(x,ty,':b',x,-ty,':b');
    hold on;
    if (BSF==0 && HPF==0)
        plot(real(DP),imag(DP),'bx','MarkerSize',12);hold on;
        plot(real(DZ),imag(DZ),'bo','MarkerSize',12);grid on;
    end
    if (BSF==1 || HPF==1)
        plot(real(DP_1),imag(DP_1),'bx','MarkerSize',12);hold on;
        plot(real(DP_2),imag(DP_2),'bo','MarkerSize',12);grid on;
    end
    title('Digital Pole-Zero Diagram');
    xlabel('Real Part');
    ylabel('Imaginary Part');
    axis square;
    
    %% Digital Magnitude and Phase Response
    
    theta = 0:0.001:pi;
    k = 0:2;
    if (BSF==0 && HPF==0)
        H = b*exp(-1i*k'*theta)./(a*exp(-1i*k'*theta));
    end
    if (BSF==1 || HPF==1)
        H = a_2*exp(-1i*k'*theta)./(a_1*exp(-1i*k'*theta));
    end
    if (d(1)==0 && d(2)==0)
        b0 = 1/H(1);
    end
    if (d(1)>0 && d(3)==0)
        b0 = abs(1/H(end));
    end
    if (d(2)>0)
        if (d(3)>0)
            b0 = sum(a)/sum(b);
        end
        if (d(3)==0)
            b0 = 1/max(abs(H));
        end
    end
    if (d(1)>0 && d(3)>0)
        b0 = 1/max(abs(H));
    end
    % b2 = b0*b;
    H = b0*H;
    figure(6);
    plot(theta/pi,abs(H));grid on;
    xlabel('\theta/\pi');ylabel('|H(\theta)|');
    title('Digital Magnitude Response (Linear)');
    figure(7);
    plot(theta/pi,20*log10(abs(H)));grid on;
    xlabel('\theta/\pi');ylabel('|H(\theta)|');
    title('Digital Magnitude Response (dB)');
    figure(8);
    plot(theta/pi,(180/pi)*angle(H));grid on;
    xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (deg)');
    title('Digital Phase Response');
    
    %% Reported Data
    
    if (d(1)==0 && d(2)==0)
        disp('H(s) is a Low-Pass Filter.');
        b0
        disp('Digital Poles:');DP
        disp('Digital Zeros:');DZ
        disp('Coeffecients of H(z) numerator:');
        b
        disp('Coeffecients of H(z) numerator multiplied by b0:');
        b0*b
        disp('Coeffecients of H(z) denominator:');
        a
    end
    
    if (d(1)>0 && d(3)==0)
        disp('H(s) is a High-Pass Filter.');
        b0
        DP = DP_1;
        DZ = DP_2;
        disp('Digital Poles:');DP
        disp('Digital Zeros:');DZ
        disp('Coeffecients of H(z) numerator:');
        a_2
        disp('Coeffecients of H(z) numerator multiplied by b0:');
        b0*a_2
        disp('Coeffecients of H(z) denominator:');
        a_1
    end
    
    if (d(2)>0)
        disp('H(s) is a Band-Pass Filter.');
        b0
        disp('Digital Poles:');DP
        disp('Digital Zeros:');DZ
        disp('Coeffecients of H(z) numerator:');
        b
        disp('Coeffecients of H(z) numerator multiplied by b0:');
        b0*b
        disp('Coeffecients of H(z) denominator:');
        a
    end
    
    if (d(1)>0 && d(3)>0)
        disp('H(s) is a Band-Stop Filter.');
        b0
        DP = DP_1;
        DZ = DP_2;
        disp('Digital Poles:');DP
        disp('Digital Zeros:');DZ
        disp('Coeffecients of H(z) numerator:');
        a_2
        disp('Coeffecients of H(z) numerator multiplied by b0:');
        b0*a_2
        disp('Coeffecients of H(z) denominator:');
        a_1
    end
    disp('First term is a constant, second term is z^-1, third term is z^-2');
    
end

if (Approx_Select == 2)
    
    LPF = 0;
    HPF = 0;
    BPF = 0;
    BSF = 0;
    
    if (num(1)==0 && num(3)>0)
        LPF = 1;
    end
    if (num(1)>0 && num(3)==0)
        HPF = 1;
    end
    if (num(2)>0)
        BPF = 1;
    end
    if (num(1)>0 && num(3)>0)
        BSF = 1;
    end
    
    if (LPF==1)
        disp('LPF');
        s = roots(den);
        n = length(s);
        sz = roots(num);
        sz = [sz;realmax*ones(1,n-length(sz))'];
        sz_actual = sz;
        index_infinity = find(sz>1000);
        %%sz(index_infinity) = 0;
        sz_actual(index_infinity) = Inf;
        z = (2*fs+sz)./(2*fs-sz);
        p = (2*fs+s)./(2*fs-s);
        b = poly(z);
        a = poly(p);
        b0 = abs(sum(a)/sum(b))*num(end)/den(end);
        b1 = b0*b;
        theta = 0:0.002:pi;
        Q = exp(-1i*(0:2)'*theta);
        H = b0*b*Q./(a*Q);
        
        %% Analog Pole-Zero Diagram
        figure(1);
        plot(real(s),imag(s),'bx','MarkerSize',12);hold on;
        plot(real(sz),imag(sz),'bo','MarkerSize',12);grid on;
        title('Analog Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        
        
        %% Magnitude and Phase Response
        figure(2);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Magnitude Response (Linear)');
        figure(3);
        plot(theta/pi,20*log10(abs(H)));grid on;
        xlabel('\theta/\pi');ylabel('|20*log10(H(\theta))|');
        title('Magnitude Response (dB)');
        figure(4);plot(theta/pi,180/pi*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (degrees)');
        title('Phase Response');
        H0 = abs(H(1));
        
        %% Digital Pole-Zero Diagram
        x = -1:0.01:1;
        ty = sqrt(1-x.^2);
        figure(5);
        plot(x,ty,':b',x,-ty,':b');
        hold on;
        plot(real(p),imag(p),'bx','MarkerSize',12);hold on;
        plot(real(z),imag(z),'bo','MarkerSize',12);grid on;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');ylabel('Imaginary Part');
        axis square;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        
        %% Digital Magnitude and Phase Response
        k = 0:2;
        figure(6);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (Linear)');
        figure(7);
        plot(theta/pi,20*log10(abs(H)));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (dB)');
        figure(8);
        plot(theta/pi,(180/pi)*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (deg)');
        title('Digital Phase Response');
        
        %% Reported Data
        b0
        disp('Analog Poles:');
        s                                                          %#ok<*NOPTS>
        disp('Analog Zeros:');
        sz_actual
        disp('Digital Poles:');
        p
        disp('Digital Zeros:');
        z
        %disp('Difference Equation:')
        %str = ['y(k)=',num2str(-a(2)),'y(k-1)-',num2str(-a(3)),'
        disp('Coeffecients of H(z) numerator:');
        b                                                     %#ok<NOPTS>
        disp('Coeffecients of H(z) numerator multiplied by b0:');
        b0*b                                                  %#ok<NOPTS>
        disp('Coeffecients of H(z) denominator:');
        a                                                     %#ok<NOPTS>
        disp('First term is a constant, second term is z^-1, third term is z^-2');
        
    end
    
    if (HPF==1)
        disp('HPF');
        s = roots(den);
        sz = roots(num);
        sz_actual = sz;
        index_infinity = find(sz>1000);
        %%sz(index_infinity) = 0;
        sz_actual(index_infinity) = Inf;
        z = (2*fs+sz)./(2*fs-sz);
        p = (2*fs+s)./(2*fs-s);
        b = poly(z);
        a = poly(p);
        L = length(a);
        n = 0:L-1;
        b0 = abs(sum((-1).^n.*a)/sum((-1).^n.*b));
        theta = 0:0.002:pi;
        Q = exp(-1i*(0:2)'*theta);
        H = b0*b*Q./(a*Q);
        
        figure(1);
        plot(real(s),imag(s),'bx','MarkerSize',12);hold on;
        plot(real(sz),imag(sz),'bo','MarkerSize',12);grid on;
        title('Analog Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        
        figure(2);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Magnitude Response (Linear)');
        Mag = 20*log10(abs(H));
        Th = -100;
        Mag = (Mag >= Th).*Mag + Th*(Mag < Th);
        figure(3);
        plot(theta/pi,Mag);grid on;
        xlabel('\theta/\pi');ylabel('|20*log10(H(\theta))/\pi');
        title('Magnitude Response (dB)');
        figure(4);
        plot(theta/pi,180/pi*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta) (deg)');
        title('Magnitude Response (dB)');
        H1 = abs(H(length(H)));
        
        %% Digital Pole-Zero Diagram
        x = -1:0.01:1;
        ty = sqrt(1-x.^2);
        figure(5);
        plot(x,ty,':b',x,-ty,':b');
        hold on;
        plot(real(p),imag(p),'bx','MarkerSize',12);hold on;
        plot(real(z),imag(z),'bo','MarkerSize',12);grid on;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');ylabel('Imaginary Part');
        axis square;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        
        %% Digital Magnitude and Phase Response
        k = 0:2;
        figure(6);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (Linear)');
        figure(7);
        plot(theta/pi,20*log10(abs(H)));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (dB)');
        figure(8);
        plot(theta/pi,(180/pi)*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (deg)');
        title('Digital Phase Response');
        
        %% Reported Data
        b0
        disp('Analog Poles:');
        s                                                          %#ok<*NOPTS>
        disp('Analog Zeros:');
        sz_actual
        disp('Digital Poles:');
        p
        disp('Digital Zeros:');
        z
        %disp('Difference Equation:')
        %str = ['y(k)=',num2str(-a(2)),'y(k-1)-',num2str(-a(3)),'
        disp('Coeffecients of H(z) numerator:');
        b                                                     %#ok<NOPTS>
        disp('Coeffecients of H(z) numerator multiplied by b0:');
        b0*b                                                  %#ok<NOPTS>
        disp('Coeffecients of H(z) denominator:');
        a                                                     %#ok<NOPTS>
        disp('First term is a constant, second term is z^-1, third term is z^-2');
        
    end
    
    if (BPF==1)
        disp('BPF');
        s = roots(den).';
        n = length(s);
        sz = roots(num).';
        sz = [sz;realmax*ones(1,n-length(sz))'];
        sz_actual = sz;
        index_infinity = find(sz>1000);
        %sz(index_infinity) = 0;
        sz_actual(index_infinity) = Inf;
        z = (2*fs+sz)./(2*fs-sz);
        p = (2*fs+s)./(2*fs-s);
        b = poly(z);
        a = poly(p);
        wp = sqrt(den(3));
        qp = 2*atan(wp/(2*fs));
        qa = exp(1i*qp*[2 1 0]);
        b0 = abs(sum(qa.*a)/sum(qa.*b));
        theta = 0:0.002:pi;
        Q = exp(-1i*(0:2)'*theta);
        H = (b0*b*Q)./(a*Q);
        n1 = 1:201;
        %x = (n1-101)/100;
        figure(1);
        plot(real(s),imag(s),'bx','MarkerSize',12);hold on;
        plot(real(sz),imag(sz),'bo','MarkerSize',12);grid on;
        title('Analog Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        figure(2);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Magnitude Response (Linear)');
        figure(3);
        plot(theta/pi,20*log10(abs(H)));grid on;
        xlabel('\theta/\pi');ylabel('|20*log10(H(\theta))|');
        title('Magnitude Response (dB)');
        figure(4);plot(theta/pi,180/pi*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (degrees)');
        title('Phase Response');
        Q1 = exp(-1i*(0:2)'*qp);
        H1 = abs(b0*b*Q1./(a*Q1));
        
        %% Digital Pole-Zero Diagram
        x = -1:0.01:1;
        ty = sqrt(1-x.^2);
        figure(5);
        plot(x,ty,':b',x,-ty,':b');
        hold on;
        plot(real(p),imag(p),'bx','MarkerSize',12);hold on;
        plot(real(z),imag(z),'bo','MarkerSize',12);grid on;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');ylabel('Imaginary Part');
        axis square;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        
        %% Digital Magnitude and Phase Response
        k = 0:2;
        figure(6);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (Linear)');
        figure(7);
        plot(theta/pi,20*log10(abs(H)));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (dB)');
        figure(8);
        plot(theta/pi,(180/pi)*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (deg)');
        title('Digital Phase Response');
        
        %% Reported Data
        b0
        disp('Analog Poles:');
        s                                                          %#ok<*NOPTS>
        disp('Analog Zeros:');
        sz_actual
        disp('Digital Poles:');
        p
        disp('Digital Zeros:');
        z
        %disp('Difference Equation:')
        %str = ['y(k)=',num2str(-a(2)),'y(k-1)-',num2str(-a(3)),'
        disp('Analog filter maximum magnitude (rad/s):');
        wp
        disp('This frequency transformed to (rad):');
        qp
        disp('Coeffecients of H(z) numerator:');
        b                                                     %#ok<NOPTS>
        disp('Coeffecients of H(z) numerator multiplied by b0:');
        b0*b                                                  %#ok<NOPTS>
        disp('Coeffecients of H(z) denominator:');
        a                                                     %#ok<NOPTS>
        disp('First term is a constant, second term is z^-1, third term is z^-2');
        
        
    end
    
    if (BSF==1)
        disp('BSF');
        s = roots(den).';
        sz = roots(num).';
        sz_actual = sz;
        index_infinity = find(sz>1000);
        %%sz(index_infinity) = 0;
        sz_actual(index_infinity) = Inf;
        z = (2*fs+sz)./(2*fs-sz);
        p = (2*fs+s)./(2*fs-s);
        b = poly(z);
        a = poly(p);
        b0 = abs((sum(a)/sum(b))*num(3)/den(3));
        theta = 0:0.002:pi;
        Q = exp(-1i*(0:2)'*theta);
        H = (b0*b*Q)./(a*Q);
        
        figure(1);
        plot(real(s),imag(s),'bx','MarkerSize',12);hold on;
        plot(real(sz),imag(sz),'bo','MarkerSize',12);grid on;
        title('Analog Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        figure(2);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Magnitude Response (Linear)');
        figure(3);
        plot(theta/pi,20*log10(abs(H)));grid on;
        xlabel('\theta/\pi');ylabel('|20*log10(H(\theta))|');
        title('Magnitude Response (dB)');
        figure(4);plot(theta/pi,180/pi*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (degrees)');
        title('Phase Response');
        H0 = H(1);
        
        %% Digital Pole-Zero Diagram
        x = -1:0.01:1;
        ty = sqrt(1-x.^2);
        figure(5);
        plot(x,ty,':b',x,-ty,':b');
        hold on;
        plot(real(p),imag(p),'bx','MarkerSize',12);hold on;
        plot(real(z),imag(z),'bo','MarkerSize',12);grid on;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');ylabel('Imaginary Part');
        axis square;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        
        %% Digital Magnitude and Phase Response
        k = 0:2;
        figure(6);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (Linear)');
        figure(7);
        plot(theta/pi,20*log10(abs(H)));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (dB)');
        figure(8);
        plot(theta/pi,(180/pi)*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (deg)');
        title('Digital Phase Response');
        
        %% Reported Data
        b0
        disp('Analog Poles:');
        s                                                          %#ok<*NOPTS>
        disp('Analog Zeros:');
        sz_actual
        disp('Digital Poles:');
        p
        disp('Digital Zeros:');
        z
        %disp('Difference Equation:')
        %str = ['y(k)=',num2str(-a(2)),'y(k-1)-',num2str(-a(3)),'
        %     disp('Analog filter maximum magnitude (rad/s):');
        %     wp
        %     disp('This frequency transformed to (rad):');
        %     qp
        disp('Coeffecients of H(z) numerator:');
        b                                                     %#ok<NOPTS>
        disp('Coeffecients of H(z) numerator multiplied by b0:');
        b0*b                                                  %#ok<NOPTS>
        disp('Coeffecients of H(z) denominator:');
        a                                                     %#ok<NOPTS>
        disp('First term is a constant, second term is z^-1, third term is z^-2');
        
    end
    
end

if (Approx_Select == 3)
    
    LPF = 0;
    HPF = 0;
    BPF = 0;
    BSF = 0;
    
    if (num(1)==0 && num(3)>0)
        LPF = 1;
    end
    if (num(1)>0 && num(3)==0)
        HPF = 1;
    end
    if (num(2)>0)
        BPF = 1;
    end
    if (num(1)>0 && num(3)>0)
        BSF = 1;
    end
    
    if (LPF==1)
        s = roots(den);
        sz = roots(num);
        sz = [sz,-realmax*ones(1,length(s)-length(sz))];
        sz_actual = sz;
        index_infinity = find(abs(sz)>1000);
        %%sz(index_infinity) = 0;
        sz_actual(index_infinity) = Inf;                           %#ok<*FNDSB>
        z = 1./(1-sz/fs);
        z_actual = z.^2;
        p = 1./(1-s/fs);
        b = fliplr(poly(z));
        a = poly(p);
        b0 = abs((sum(a)/sum(b))*num(3)/den(3));
        theta = 0:0.01:pi;
        Q = exp(-1i*(0:2)'*theta);
        H = b0*b*Q./(a*Q);
        
        %% Analog Pole-Zero Diagram
        figure(1);
        plot(real(s),imag(s),'bx','MarkerSize',12);hold on;
        plot(real(sz),imag(sz),'bo','MarkerSize',12);grid on;
        title('Analog Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        
        
        %% Magnitude and Phase Response
        figure(2);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Magnitude Response (Linear)');
        figure(3);
        plot(theta/pi,20*log10(abs(H)));grid on;
        xlabel('\theta/\pi');ylabel('|20*log10(H(\theta))|');
        title('Magnitude Response (dB)');
        figure(4);plot(theta/pi,180/pi*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (degrees)');
        title('Phase Response');
        H0 = abs(H(1));
        
        %% Digital Pole-Zero Diagram
        x = -1:0.01:1;
        ty = sqrt(1-x.^2);
        figure(5);
        plot(x,ty,':b',x,-ty,':b');
        hold on;
        plot(real(p),imag(p),'bx','MarkerSize',12);hold on;
        plot(real(z),imag(z),'bo','MarkerSize',12);grid on;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');ylabel('Imaginary Part');
        axis square;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        
        %% Digital Magnitude and Phase Response
        k = 0:2;
        figure(6);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (Linear)');
        figure(7);
        plot(theta/pi,20*log10(abs(H)));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (dB)');
        figure(8);
        plot(theta/pi,(180/pi)*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (deg)');
        title('Digital Phase Response');
        
        %% Reported Data
        b0
        disp('Analog Poles:');
        s                                                          %#ok<*NOPTS>
        disp('Analog Zeros:');
        sz_actual
        disp('Digital Poles:');
        p
        disp('Digital Zeros:');
        z_actual
        %disp('Difference Equation:')
        %str = ['y(k)=',num2str(-a(2)),'y(k-1)-',num2str(-a(3)),'
        disp('Coeffecients of H(z) numerator:');
        b                                                     %#ok<NOPTS>
        disp('Coeffecients of H(z) numerator multiplied by b0:');
        b0*b                                                  %#ok<NOPTS>
        disp('Coeffecients of H(z) denominator:');
        a                                                     %#ok<NOPTS>
        disp('First term is a constant, second term is z^-1, third term is z^-2');
        
    end
    
    if (HPF==1)
        sz = roots(num);
        s = roots(den);
        sz_actual = sz;
        index_infinity = find(abs(sz)>1000);
        %%sz(index_infinity) = 0;
        sz_actual(index_infinity) = Inf;                           %#ok<*FNDSB>
        z = 1./(1-sz/fs);
        z_actual = z.^2;
        p = 1./(1-s/fs);
        b = poly(z);
        a = poly(p);
        L = length(a);
        n = 0:L-1;
        b0 = abs(sum((-1).^n.*a)/sum((-1).^n.*b));
        theta = 0:0.002:pi;
        Q = exp(-1i*(0:2)'*theta);
        H = b0*b*Q./(a*Q);
        
        %% Analog Pole-Zero Diagram
        figure(1);
        plot(real(s),imag(s),'bx','MarkerSize',12);hold on;
        plot(real(sz),imag(sz),'bo','MarkerSize',12);grid on;
        title('Analog Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        
        
        %% Magnitude and Phase Response
        figure(2);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Magnitude Response (Linear)');
        figure(3);
        plot(theta/pi,20*log10(abs(H)));grid on;
        xlabel('\theta/\pi');ylabel('|20*log10(H(\theta))|');
        title('Magnitude Response (dB)');
        figure(4);plot(theta/pi,180/pi*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (degrees)');
        title('Phase Response');
        H1 = abs(H(length(H)));
        
        %% Digital Pole-Zero Diagram
        x = -1:0.01:1;
        ty = sqrt(1-x.^2);
        figure(5);
        plot(x,ty,':b',x,-ty,':b');
        hold on;
        plot(real(p),imag(p),'bx','MarkerSize',12);hold on;
        plot(real(z),imag(z),'bo','MarkerSize',12);grid on;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');ylabel('Imaginary Part');
        axis square;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        
        %% Digital Magnitude and Phase Response
        k = 0:2;
        figure(6);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (Linear)');
        figure(7);
        plot(theta/pi,20*log10(abs(H)));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (dB)');
        figure(8);
        plot(theta/pi,(180/pi)*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (deg)');
        title('Digital Phase Response');
        
        %% Reported Data
        b0
        disp('Analog Poles:');
        s                                                          %#ok<*NOPTS>
        disp('Analog Zeros:');
        sz_actual
        disp('Digital Poles:');
        p
        disp('Digital Zeros:');
        z_actual
        %disp('Difference Equation:')
        %str = ['y(k)=',num2str(-a(2)),'y(k-1)-',num2str(-a(3)),'
        disp('Coeffecients of H(z) numerator:');
        b                                                     %#ok<NOPTS>
        disp('Coeffecients of H(z) numerator multiplied by b0:');
        b0*b                                                  %#ok<NOPTS>
        disp('Coeffecients of H(z) denominator:');
        a                                                     %#ok<NOPTS>
        disp('First term is a constant, second term is z^-1, third term is z^-2');
    end
    
    if (BPF==1)
        sz = roots(num);
        if (isempty(sz) )
            sz=[-inf -inf];
        end
        if (length(sz)==1)
            sz=[sz -inf]; end
        sz;
        s = roots(den);
        sz_actual = sz;
        index_infinity = find(abs(sz)>1000);
        %%sz(index_infinity) = 0;
        sz_actual(index_infinity) = Inf;                           %#ok<*FNDSB>
        z = 1./(1-sz/fs);
        z_actual = z.^2;
        p = 1./(1-s/fs);
        b = poly(z);
        a = poly(p);
        k = 0:2;
        a1 = sum(p);
        b1 = prod(p);
        c1 = z(1);
        thetap = acos((4*b1+4*b1*c1^2-4*sqrt(b1^2+b1^2*c1^4+b1*c1^2+a1^2*b1*c1^2+b1^3*c1^2-a1*b1*c1-a1*b1^2*c1-a1*b1*c1^3-a1*b1^2*c1^3))/(8*b1*c1));
        H2 = b*exp(-1i*k'*thetap)./(a*exp(-1i*k'*thetap));
        b0 = 1/abs(H2);
        theta = 0:0.001:pi;
        Q = exp(-1i*k'*theta);
        H = b0*b*Q./(a*Q);
        
        figure(1);
        plot(real(s),imag(s),'bx','MarkerSize',12);hold on;
        plot(real(sz),imag(sz),'bo','MarkerSize',12);grid on;
        title('Analog Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        
        
        figure(2);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Magnitude Response (Linear)');
        Mag = 20*log10(abs(H));
        Th = -100;
        Mag = (Mag>=Th).*Mag+Th*(Mag<Th);
        figure(3);
        plot(theta/pi,Mag);grid on;
        xlabel('\theta/\pi');ylabel('|20*log10(H(\theta))|');
        title('Magnitude Response (dB)');
        figure(4);plot(theta/pi,180/pi*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (degrees)');
        title('Phase Response');
        
        %% Digital Pole-Zero Diagram
        x = -1:0.01:1;
        ty = sqrt(1-x.^2);
        figure(5);
        plot(x,ty,':b',x,-ty,':b');
        hold on;
        plot(real(p),imag(p),'bx','MarkerSize',12);hold on;
        plot(real(z),imag(z),'bo','MarkerSize',12);grid on;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');ylabel('Imaginary Part');
        axis square;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        
        %% Digital Magnitude and Phase Response
        k = 0:2;
        figure(6);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (Linear)');
        figure(7);
        plot(theta/pi,20*log10(abs(H)));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (dB)');
        figure(8);
        plot(theta/pi,(180/pi)*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (deg)');
        title('Digital Phase Response');
        
        %% Reported Data
        b0
        disp('Analog Poles:');
        s                                                          %#ok<*NOPTS>
        disp('Analog Zeros:');
        sz_actual
        disp('Digital Poles:');
        p
        disp('Digital Zeros:');
        z
        %disp('Difference Equation:')
        %str = ['y(k)=',num2str(-a(2)),'y(k-1)-',num2str(-a(3)),'
        %     disp('Analog filter maximum magnitude (rad/s):');
        %     wp
        %     disp('This frequency transformed to (rad):');
        %     qp
        disp('Coeffecients of H(z) numerator:');
        b                                                     %#ok<NOPTS>
        disp('Coeffecients of H(z) numerator multiplied by b0:');
        b0*b                                                  %#ok<NOPTS>
        disp('Coeffecients of H(z) denominator:');
        a                                                     %#ok<NOPTS>
        disp('First term is a constant, second term is z^-1, third term is z^-2');
    end
    
    if (BSF==1)
        sz = roots(num);
        s = roots(den);
        sz_actual = sz;
        index_infinity = find(abs(sz)>1000);
        %%sz(index_infinity) = 0;
        sz_actual(index_infinity) = Inf;                           %#ok<*FNDSB>
        z = 1./(1-sz/fs);
        z_actual = z.^2;
        p = 1./(1-s/fs);
        b = poly(z);
        a = poly(p);
        b0 = abs((sum(a)/sum(b))*num(3)/den(3));
        theta = 0:0.002:pi;
        Q = exp(-1i*(0:2)'*theta);
        H= b0*b*Q./(a*Q);
        
        figure(1);
        plot(real(s),imag(s),'bx','MarkerSize',12);hold on;
        plot(real(sz),imag(sz),'bo','MarkerSize',12);grid on;
        title('Analog Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        figure(2);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Magnitude Response (Linear)');
        figure(3);
        plot(theta/pi,20*log10(abs(H)));grid on;
        xlabel('\theta/\pi');ylabel('|20*log10(H(\theta))|');
        title('Magnitude Response (dB)');
        figure(4);plot(theta/pi,180/pi*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (degrees)');
        title('Phase Response');
        H0 = H(1);
        
        %% Digital Pole-Zero Diagram
        x = -1:0.01:1;
        ty = sqrt(1-x.^2);
        figure(5);
        plot(x,ty,':b',x,-ty,':b');
        hold on;
        plot(real(p),imag(p),'bx','MarkerSize',12);hold on;
        plot(real(z),imag(z),'bo','MarkerSize',12);grid on;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');ylabel('Imaginary Part');
        axis square;
        title('Digital Pole-Zero Diagram');
        xlabel('Real Part');
        ylabel('Imaginary Part');
        axis square;
        
        %% Digital Magnitude and Phase Response
        k = 0:2;
        figure(6);
        plot(theta/pi,abs(H));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (Linear)');
        figure(7);
        plot(theta/pi,20*log10(abs(H)));grid on;
        xlabel('\theta/\pi');ylabel('|H(\theta)|');
        title('Digital Magnitude Response (dB)');
        figure(8);
        plot(theta/pi,(180/pi)*angle(H));grid on;
        xlabel('\theta/\pi');ylabel('Angle(H(\theta)) (deg)');
        title('Digital Phase Response');
        
        %% Reported Data
        b0
        disp('Analog Poles:');
        s                                                          %#ok<*NOPTS>
        disp('Analog Zeros:');
        sz_actual
        disp('Digital Poles:');
        p
        disp('Digital Zeros:');
        z
        %disp('Difference Equation:')
        %str = ['y(k)=',num2str(-a(2)),'y(k-1)-',num2str(-a(3)),'
        %     disp('Analog filter maximum magnitude (rad/s):');
        %     wp
        %     disp('This frequency transformed to (rad):');
        %     qp
        disp('Coeffecients of H(z) numerator:');
        b                                                     %#ok<NOPTS>
        disp('Coeffecients of H(z) numerator multiplied by b0:');
        b0*b                                                  %#ok<NOPTS>
        disp('Coeffecients of H(z) denominator:');
        a                                                     %#ok<NOPTS>
        disp('First term is a constant, second term is z^-1, third term is z^-2');
        
    end
    
end
