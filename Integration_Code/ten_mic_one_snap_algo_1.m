function [mj, si, ti, eps_best, theta_best] = ten_mic_one_snap_algo_1(tij, M, err, thresh, peaks)
    theta_orig = 340 * tij;
    M_s = 340 * M;
    M_trans = M_s';
    num_peaks = size(M,2);
    num_mic = size(M,1);
    theta_est = zeros(size(M_trans));
    eps_best = Inf;
    dim = 3;
    
    for i1 = 1:num_peaks
        theta_est(i1, 1) = M_trans(i1, 1);
        for i2 = 1:num_peaks
            if abs(M_trans(i2, 2) - theta_est(i1,1)) < err
                theta_est(i2, 2) = M_trans(i2, 2);
            else
                continue
            end
            for i3 = 1:num_peaks
                if abs(M_trans(i3, 3) - theta_est(i1,1)) < err
                    theta_est(i3, 3) = M_trans(i3, 3);
                else
                    continue
                end
                for i4 = 1:num_peaks
                    if abs(M_trans(i4, 4) - theta_est(i1,1)) < err
                        theta_est(i4, 4) = M_trans(i4, 4);
                    else
                        continue
                    end
                    for i5 = 1:num_peaks
                        if abs(M_trans(i5, 5) - theta_est(i1,1)) < err
                            theta_est(i5, 5) = M_trans(i5, 5);
                        else
                            continue
                        end
                        for i6 = 1:num_peaks
                            if abs(M_trans(i6, 6) - theta_est(i1,1)) < err
                                theta_est(i6, 6) = M_trans(i6, 6);
                            else
                                continue
                            end
                            for i7 = 1:num_peaks
                                if abs(M_trans(i7, 7) - theta_est(i1,1)) < err
                                    theta_est(i7, 7) = M_trans(i7, 7);
                                else
                                    continue
                                end
                                for i8 = 1:num_peaks
                                    if abs(M_trans(i8, 8) - theta_est(i1,1)) < err
                                        theta_est(i8, 8) = M_trans(i8, 8);
                                    else
                                        continue
                                    end
                                    for i9 = 1:num_peaks
                                        if abs(M_trans(i9, 9) - theta_est(i1,1)) < err
                                            theta_est(i9, 9) = M_trans(i9, 9);
                                        else
                                            continue
                                        end
                                        for i10 = 1:num_peaks
                                            if abs(M_trans(i10, 10) - theta_est(i1,1)) < err
                                                theta_est(i10, 10) = M_trans(i10, 10);
                                            else
                                                continue
                                            end
                                            
                                            for ii = 1:7
                                                [~, echoScores] = sort_echoes(D, peaks(:, ii), dim, num_mic);
                                                for jj = 1:len(echoScores)
                                                    if echoScores(jj) > thresh
                                                        break
                                                    end
                                                end
                                            end
                                            
                                            [ti] = calctod_2D_sumanth(theta_est);
                                            [Si_cap, Mj_cap] = refine_absol_2D_sumanth_tdoa(theta_est,ti);

                                            R = [1 0; 0 1];
                                            t = [0 0];

                                            [Si,Mj] = compute_locations_2D_sumanth(Si_cap,Mj_cap,R,t);
                                            [si,mj] = obtain_source_mic_locations(Si,Mj);

                                            eps = norm(theta_orig - theta_est,'fro');

                                            if eps < eps_best
                                                eps_best = eps;
                                                theta_best = theta_est;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end