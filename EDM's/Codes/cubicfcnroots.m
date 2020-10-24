function [x,nroot]=cubicfcnroots(a, b, c, d)%For solving equations of the form a*x^3 + b*x^2 + c*x + d=0

      if (a == 0.)%If a is 0, it remains a quadratic equation
        [x, nroot] = quadfcn(b, c, d);
        return
      end

      nroot = 3;
      p  = c/a - b*b/a/a/3. ;%Calculate p and q
      q  = (2.*b*b*b/a/a/a - 9.*b*c/a/a + 27.*d/a) / 27. ;

      D = p*p*p/27. + q*q/4. ;%Calculate the Discriminant
      
      %If D<0, three real roots
      %If D=0, three real roots, of which two are equal
      %If D>0, one real and two imaginary roots
      
     % For D>0 and D=0,
%           Calculate u and v
%           u = cubic_root(-q/2 + sqrt(D))
%           v = cubic_root(-q/2 - sqrt(D))
%           Find the three transformed roots
%           y1 = u + v
%           y2 = -(u+v)/2 + i (u-v)*sqrt(3)/2
%           y3 = -(u+v)/2 - i (u-v)*sqrt(3)/2
%   Step 4  Finally, find the three roots
%           x = y - b/a/3

    %For D<0, a trigonometric formulation is more convenient
    %y1 =  2 * sqrt(|p|/3) * cos(phi/3) 
    %y2 = -2 * sqrt(|p|/3) * cos((phi+pi)/3)
    %y3 = -2 * sqrt(|p|/3) * cos((phi-pi)/3)
    % where phi = acos(-q/2/sqrt(|p|**3/27)) & pi  = 3.141592654...
      if (D < 0.)
        phi = acos(-q/2./sqrt(abs(p*p*p)/27.));
        temp1=2.*sqrt(abs(p)/3.);
        y1 =  temp1*cos(phi/3.);
        y2 = -temp1*cos((phi+pi)/3.);
        y3 = -temp1*cos((phi-pi)/3.);
      else
        temp1 = -q/2. + sqrt(D);  % For D>0 and D=0,
        temp2 = -q/2. - sqrt(D);%Calculate u and v
        u = abs(temp1)^(1./3.);%u = cubic_root(-q/2 + sqrt(D))
        v = abs(temp2)^(1./3.); %v = cubic_root(-q/2 - sqrt(D))
        if (temp1 < 0.) %Find the three transformed roots
            u=-u; %y1 = u + v
        end
        if (temp2 < 0.) %y2 = -(u+v)/2 + i (u-v)*sqrt(3)/2 
            v=-v; %y3 = -(u+v)/2 - i (u-v)*sqrt(3)/2
        end
        y1  = u + v;
        y2r = -(u+v)/2.;
        y2i =  (u-v)*sqrt(3.)/2.;
      end

      temp1 = b/a/3.;
      y1 = y1-temp1;
      if (D < 0.)
        y2 = y2-temp1;
        y3 = y3-temp1;
      else
        y2r=y2r-temp1;
      end

      if (D < 0.)
        x(1) = y1;
        x(2) = y2;
        x(3) = y3;
      elseif (D == 0.)
        x(1) =  y1;
        x(2) = y2r;
        x(3) = y2r;
      else
        x(1) = y1;
        x(2) = y2r + y2i*1i;
        x(3) = y2r - y2i*1i;
      end
      x = x(:);
end