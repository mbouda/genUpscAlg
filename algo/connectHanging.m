function [hookEqs,nHooks]=connectHanging(hangLinks,prob,parents)

   keyboard     
        
        nHL=size(hangLinks,1);
        
        isConn=false(nHL,1);
        tops=zeros(nHL,1);
        hits=zeros(nHL,1);
        pars=parents(hangLinks);
        
        while any(~isConn)
            
            [~,imp]=max(pars);
            
            [touches,hit]=srchTarg(pars(imp),parents,cat(1,prob.targ,hangLinks(~imp)),cat(1,prob.bots,prob.terms,hangLinks(imp)));
            if touches
                isConn(imp)=true;
                hits(imp)=hit;
                tops(imp)=pars(imp);
                pars(imp)=0;
            else
                pars(imp)=parents(pars(imp));
            end
        end
        
        for i=1:nHL
            %need to know what it connects to...
            %merge connecteds: 
                %give them single parent, 
                %make both hangers ineligible to end the search
                
            %will then need to continue search until all hangers are
            %connected to an actual target, though keep track of the
            %individual equations on the way...
                %likely a while loop around the one we have above
                        
        end
        
        %also: do we need to advance by layer (likely not: not adding psiBarX equations)? 
        %should count d.o.f.? as in, extra layers... the tops... 
        %can define correct target number of d.o.f.s? Do have info on
        %bots/targs & nLayers in here... maybe need # extraEqs as well...
        
        
end