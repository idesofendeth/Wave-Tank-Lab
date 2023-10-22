% %% Run both the horizontal and vertical case and save the concatenated data
for i = 1:2
    CaseNr = i;
    
    switch CaseNr
        case 1
%             horiz_vertic_sametime;
%             R=R_members;
%Horizontal_Radii;
% radii_ismember_mixed_9Aug;
% calc_dR_and_Lambda_mixed_9Aug
Horizontal_Radii;
            calc_dR_and_Lambda_test25Aug;
            save celldata.mat lam_i c_exp_lead c_exp_trail
            %save DATA22_horizontal_Marble_6cm.mat EXPDATA
            save DATA22_horizontal_Droplet_6cm.mat EXPDATA
            x_horizontal_rhs = x_crests_rhs; %Save x-data with new name to be able to animate horizontally and vertically together
            x_horizontal_lhs = x_crests_lhs;
        case 2
%             horiz_vertic_sametime;
%             R=R_members;
 Vertical_Radii;
            calc_dR_and_Lambda_test25Aug;
            %calc_dR_and_Lambda_mixed_9Aug
            %save DATA22_vertical_Marble_6cm.mat EXPDATA
            save DATA22_vertical_Droplet_6cm.mat EXPDATA
            save celldata2.mat lam_i c_exp_lead c_exp_trail
%             y_vertical_rhs = y_crests_above;
%             y_vertical_lhs = y_crests_below;
            y_vertical_rhs = x_crests_rhs;
            y_vertical_lhs = x_crests_lhs;
    end
end