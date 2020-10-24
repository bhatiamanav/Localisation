function [x,nroot]=quadfcn(a, b, c)%For solving equations of the form a*x^2 + b*x + c = 0
      if (a == 0)
        if (b == 0)
          nroot = 0;%If a and b are zero,it is a non-equation
          x = [];
        else
          nroot = 1;%else if it is a linear equation, then b*x + c = 0
          x = -c/b;
        end
      else
        nroot = 2;%Using the ShreeDharacharya formula
        D = b*b-4*a*c;
        x(1) = (-b+sqrt(D))/2/a;
       x(2) = (-b-sqrt(D))/2/a;
      end
end