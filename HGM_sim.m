function results = HGM_sim(sim)
    
    % Run simulations.
    %
    % USAGE: results = HGM_sim(sim)
    
    opts = set_opts;
    set(0, 'DefaultAxesFontName', 'Palatino');
    set(0, 'DefaultTextFontName', 'Palatino');
    
    switch sim
        
        case 'dispersion'
            
            % Weber & Crocker (1983); Hewstone & Hamberger (2000)
            
            N = 6;
            X = ones(N+2,2);
            parents = [1 1];
            Z = zeros(N+2,2);
            Z(1:N,1) = 1;
            z = {[1 0; 1 0] [1 0; 0 1] [0 1; 1 0] [0 1; 0 1]};
            x = {[1 0; 0 1] [1 1; 0 0]};
            
            for i = 1:2
                X(N+1:end,:) = x{i};
                for j = 1:length(z)
                    Z(N+1:end,:) = z{j};
                    score(j) = HGM_score(X,Z,parents,opts);
                    mu_rel(j,:,:) = HGM_infer(X,Z,opts);
                end
                results(1).p(i,:) = exp(score - logsumexp(score,2));
                for d = 1:2
                    results(1).mu_rel(i,:,d) = results(1).p(i,:)*squeeze(mu_rel(:,:,d));
                end
            end
            
            mu_rel0 = HGM_infer(X(1:N,:),Z(1:N,:),opts);
            results(1).change = results(1).mu_rel - mu_rel0';
            results(1).change = squeeze(mean(results(1).change(:,1,:),3));
            
            X(1:N,:) = repmat(linspace(0,1,N)',1,2);
            
            for i = 1:2
                X(N+1:end,:) = x{i};
                for j = 1:length(z)
                    Z(N+1:end,:) = z{j};
                    score(j) = HGM_score(X,Z,parents,opts);
                    mu_rel(j,:,:) = HGM_infer(X,Z,opts);
                end
                results(2).p(i,:) = exp(score - logsumexp(score,2));
                for d = 1:2
                    results(2).mu_rel(i,:,d) = results(2).p(i,:)*squeeze(mu_rel(:,:,d));
                end
            end
            
            mu_rel0 = HGM_infer(X(1:N,:),Z(1:N,:),opts);
            results(2).change = results(2).mu_rel - mu_rel0';
            results(2).change = squeeze(mean(results(2).change(:,1,:),3));
            
            figure;
            subplot(1,3,1);
            p1 = results(1).p*[1 1 0 0; 1 0 1 0]';
            p2 = results(2).p*[1 1 0 0; 1 0 1 0]';
            P = [p1(:,2) p2(:,2)];
            bar(P');
            set(gca,'FontSize',25,'XTickLabel',{'Low' 'High'},'XLim',[0.5 2.5],'YLim',[0 1]);
            xlabel('Variability','FontSize',25);
            ylabel('Probability','FontSize',25);
            legend({'Dispersed' 'Concentrated'},'FontSize',22,'Box','off','Location','NorthWest');
            mytitle('A','Left','FontSize',25,'FontWeight','Bold');
            
            subplot(1,3,2);
            C = [results(1).change results(2).change];
            bar(abs(C'));
            set(gca,'FontSize',25,'XTickLabel',{'Low' 'High'},'XLim',[0.5 2.5],'YLim',[-0.01 0.07]);
            xlabel('Variability','FontSize',25);
            ylabel('Absolute stereotype change','FontSize',25);
            mytitle('B','Left','FontSize',25,'FontWeight','Bold');
            
            N = 30;
            X = ones(N+2,2);
            parents = [1 1];
            Z = zeros(N+2,2);
            Z(1:N,1) = 1;
            z = {[1 0; 1 0] [1 0; 0 1] [0 1; 1 0] [0 1; 0 1]};
            x = {[1 0; 0 1] [1 1; 0 0]};
            
            for i = 1:2
                X(N+1:end,:) = x{i};
                for j = 1:length(z)
                    Z(N+1:end,:) = z{j};
                    score(j) = HGM_score(X,Z,parents,opts);
                    mu_rel(j,:,:) = HGM_infer(X,Z,opts);
                end
                results(3).p(i,:) = exp(score - logsumexp(score,2));
                for d = 1:2
                    results(3).mu_rel(i,:,d) = results(3).p(i,:)*squeeze(mu_rel(:,:,d));
                end
            end
            
            mu_rel0 = HGM_infer(X(1:N,:),Z(1:N,:),opts);
            results(3).change = results(3).mu_rel - mu_rel0';
            results(3).change = squeeze(mean(results(3).change(:,1,:),3));
            
            subplot(1,3,3);
            C = [results(1).change results(3).change];
            plot(abs(C'),'-o','LineWidth',5,'MarkerFaceColor','w','MarkerSize',12);
            set(gca,'FontSize',25,'XTick',[1 2],'XTickLabel',{'6' '30'},'XLim',[0.5 2.5],'YLim',[-0.01 0.07]);
            xlabel('Sample size','FontSize',25);
            ylabel('Absolute stereotype change','FontSize',25);
            mytitle('C','Left','FontSize',25,'FontWeight','Bold');
            
            set(gcf,'Position',[200 200 1900 450]);
            
        case 'typicality'
            
            % Kunda & Oleson (1995)
            
            N = 6;
            X{1} = [ones(N,1); 0.5];
            X{2} = [repmat([1 0],N/2,1); repmat([1 1],N/2,1); 0.5 1];
            X{3} = [repmat([1 1],N/2,1); repmat([1 1],N/2,1); 0.5 1];
            parents = [1 1];
            Z = zeros(N+1,2);
            Z(1:N,1) = 1;
            z = [1 0; 0 1];
            
            for i = 1:3
                for j = 1:size(z,1)
                    Z(N+1,:) = z(j,:);
                    score(j) = HGM_score(X{i},Z,parents,opts);
                    mu_rel(j,:,:) = HGM_infer(X{i},Z,opts);
                end
                results.p(i,:) = exp(score - logsumexp(score,2));
                results.mu_rel(i,:) = results.p(i,:)*squeeze(mu_rel(:,:,1));
                
                mu_rel0 = HGM_infer(X{i}(1:N,:),Z(1:N,:),opts);
                results.change(i,:) = results.mu_rel(i,:) - mu_rel0(:,1)';
                
                clear score mu_rel
            end
            
            subplot(1,2,1);
            p = results.p;
            bar(results.p(:,1));
            set(gca,'FontSize',25,'XTickLabel',{'None' 'Neutral' 'Typical'},'XLim',[0.5 3.5],'YLim',[0 1]);
            ylabel('Probability','FontSize',25);
            xlabel('Feature type','FontSize',25);
            mytitle('A','Left','FontSize',25,'FontWeight','Bold');
            subplot(1,2,2);
            bar(abs(results.change(:,1)));
            set(gca,'FontSize',25,'XTickLabel',{'None' 'Neutral' 'Typical'},'XLim',[0.5 3.5]);
            ylabel('Absolute stereotype change','FontSize',25);
            set(gcf,'Position',[200 200 950 450]);
            xlabel('Feature type','FontSize',25);
            mytitle('B','Left','FontSize',25,'FontWeight','Bold');
            
        case 'deviance'
            
            % Kunda & Oleson (1997)
            
            N = 6;
            X = ones(N+1,1);
            parents = [1 1];
            Z = [ones(N+1,1) zeros(N+1,1)];
            
            x = linspace(-1,1,50);
            z = [1 0; 0 1];
            
            for i = 1:length(x)
                X(N+1,1) = x(i);
                for j = 1:size(z,1)
                    Z(N+1,:) = z(j,:);
                    score(j) = HGM_score(X,Z,parents,opts);
                    mu_rel(j,:) = HGM_infer(X,Z,opts);
                end
                results.p(i,:) = exp(score - logsumexp(score,2));
                results.mu_rel(i,:) = results.p(i,:)*mu_rel;
            end
            
            figure;
            subplot(1,2,1);
            dist = repmat(1 - x',1,size(z,1));
            plot(dist,results.p(:,1),'LineWidth',5,'Color','b');
            set(gca,'FontSize',25,'YLim',[0 1]);
            ylabel('Probability','FontSize',25);
            xlabel('Degree of deviance','FontSize',25);
            mytitle('A','Left','FontSize',25,'FontWeight','Bold');
            
            subplot(1,2,2);
            mu_rel0 = HGM_infer(X(1:N,:),Z(1:N,:),opts);
            results.change = results.mu_rel - mu_rel0';
            plot(dist(:,1),abs(results.change(:,1)),'LineWidth',5,'Color','b');
            set(gca,'FontSize',25);
            xlabel('Degree of deviance','FontSize',25);
            ylabel('Absolute stereotype change','FontSize',25);
            mytitle('B','Left','FontSize',25,'FontWeight','Bold');
            
            set(gcf,'Position',[200 200 950 450]);
            
        case 'load'
            
             % Yzerbyrt et al. (1999)
            
            N = 6;
            X = [ones(N,1); 0.5];
            parents = [1 1];
            Z = [ones(N+1,1) zeros(N+1,1)];
            z = [1 0; 0 1];
            alpha = linspace(0.1,10,20);
            
            for i = 1:length(alpha)
                opts.alpha = alpha(i);
                for j = 1:size(z,1)
                    Z(N+1,:) = z(j,:);
                    score(j) = HGM_score(X,Z,parents,opts);
                    mu_rel(j,:) = HGM_infer(X,Z,opts);
                end
                results.p(i,:) = exp(score - logsumexp(score,2));
                results.mu_rel(i,:) = results.p(i,:)*mu_rel;
            end
            
            figure;
            subplot(1,2,1);
            plot(alpha,results.p(:,1),'LineWidth',5,'Color','b');
            set(gca,'FontSize',25,'YLim',[0 1]);
            ylabel('Probability','FontSize',25);
            xlabel('Concentration parameter (\alpha)','FontSize',25);
            mytitle('A','Left','FontSize',25,'FontWeight','Bold');
            
            subplot(1,2,2);
            mu_rel0 = HGM_infer(X(1:N,:),Z(1:N,:),opts);
            results.change = results.mu_rel - mu_rel0';
            plot(alpha,abs(results.change(:,1)),'LineWidth',5,'Color','b');
            set(gca,'FontSize',25);
            xlabel('Concentration parameter (\alpha)','FontSize',25);
            ylabel('Absolute stereotype change','FontSize',25);
            mytitle('B','Left','FontSize',25,'FontWeight','Bold');
            
            set(gcf,'Position',[200 200 950 450]);
            
        case 'dispersion_subgroup'
            
            N = 10;
            X = ones(N+2,2);
            parents = [1 1 2];
            Z = zeros(N+2,3);
            Z(1:N,1) = 1;
            z = {[1 0 0; 1 0 0] [1 0 0; 0 1 0] [1 0 0; 1 0 1] [0 1 0; 1 0 0] [0 1 0; 0 1 0] [0 1 0; 1 0 1]};
            x = {[1 0; 0 1] [1 1; 0 0]};
            
            for m = 1:2
                if m==2
                    Z = zeros(N+2,3);
                    Z(1:N/2,1) = 1;
                    Z(N/2+1:N,3) = 1;
                end
                for i = 1:2
                    X(N+1:end,:) = x{i};
                    for j = 1:length(z)
                        Z(N+1:end,:) = z{j};
                        score(j) = HGM_score(X,Z,parents,opts);
                        mu_rel(j,:,:) = HGM_infer(X,Z,opts);
                    end
                    results(m).p(i,:) = exp(score - logsumexp(score,2));
                    for d = 1:2
                        results(m).mu_rel(i,:,d) = results(m).p(i,:)*squeeze(mu_rel(:,:,d));
                    end
                end
                mu_rel0 = HGM_infer(X(1:N,:),Z(1:N,:),opts);
                change = results(m).mu_rel - mu_rel0';
                %results(m).change = squeeze(mean(change(:,1,:),3));
                results(m).change = squeeze(mean(change(:,1,:)+change(:,3,:),3));
            end
            
            figure;
            subplot(1,2,1);
            P(1,:) = results(1).p*[1 0 1 1 0 1]';
            P(2,:) = results(2).p*[1 0 1 1 0 1]';
            bar(P);
            set(gca,'FontSize',25,'XTickLabel',{'Subgroup' 'No subgroup'},'XLim',[0.5 2.5],'YLim',[0 1]);
            ylabel('Probability','FontSize',25);
            legend({'Dispersed' 'Concentrated'},'FontSize',22,'Box','off','Location','NorthEast');
            mytitle('A','Left','FontSize',25,'FontWeight','Bold');
            
            subplot(1,2,2);
            C = [results(1).change results(2).change];
            bar(abs(C'));
            set(gca,'FontSize',25,'XTickLabel',{'Subgroup' 'No subgroup'},'XLim',[0.5 2.5]);
            ylabel('Absolute stereotype change','FontSize',25);
            mytitle('B','Left','FontSize',25,'FontWeight','Bold');
            
            set(gcf,'Position',[200 200 950 450]);
            
        case 'ingroup'
            
            N = 10:30;
            parents = [1 2 2];
            
            for i = 1:length(N)
                X = [ones(N(i),1); zeros(N(i),1)];
                Z1 = [repmat([1 1 0],N(i),1); repmat([1 0 1],N(i),1)];
                Z2 = repmat([1 0 0],N(i)*2,1);
                score(1) = HGM_score(X,Z1,parents,opts);
                score(2) = HGM_score(X,Z2,parents,opts);
                results.p(i,:) = exp(score - logsumexp(score,2));
            end
            
            plot(N*2,results.p(:,1),'LineWidth',5,'Color','b');
            set(gca,'FontSize',25,'YLim',[0 1]);
            ylabel('Probability','FontSize',25);
            xlabel('Sample size','FontSize',25);
            
    end