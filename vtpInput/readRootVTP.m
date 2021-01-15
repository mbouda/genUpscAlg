function vtpRoot=readRootVTP(fileName)
    rootVTP = xmlread(fileName) ;


    vtkFile=rootVTP.getFirstChild;

    node=vtkFile.getFirstChild;
    polyData=node.getNextSibling;
    node=polyData.getFirstChild;
    piece=node.getNextSibling;
    pieceAttributes=piece.getAttributes;
    pieceDimensions=cell(1,2);
    for i=1:2
        pieceDimensions{i}=pieceAttributes.item(i-1);
    end
    %currently written as text in cells, number of lines, then number of points

    childNodes = piece.getChildNodes;
    numChildNodes = childNodes.getLength;
    vtpRoot=struct([]);
    for i=1:numChildNodes
        thisChild=childNodes.item(i-1);
        if thisChild.hasChildNodes
            nodeName=thisChild.getNodeName;

            eval(sprintf('vtpRoot(1).%s=struct();',nodeName))

            dataNodes = thisChild.getChildNodes;
            numDataNodes = dataNodes.getLength;
            for j=1:numDataNodes
                thisDataNode=dataNodes.item(j-1);
                if thisDataNode.hasChildNodes
                    dataAtt=thisDataNode.getAttributes;
                    nAtt=dataAtt.getLength;
                    dataAttributes=cell(nAtt,2);
                    for k=1:nAtt
                        att=dataAtt.item(k-1);
                        dataAttributes{k,1}=att.getName.toCharArray';
                        dataAttributes{k,2}=att.getValue.toCharArray';
                    end
                    data=thisDataNode.getFirstChild;
                    content=char(data.getData);
                    nums=regexp(content,'\s+','split');
                    numArray=cellfun(@(x)str2double(x),nums);
                    numArray(isnan(numArray))=[];
                    eval(sprintf('vtpRoot(1).%s(1).%s=struct();',nodeName,dataAttributes{1,2}))


                    %write to structure
                    for k=2:nAtt
                        eval(sprintf('vtpRoot(1).%s(1).%s.%s=dataAttributes{k,2};',nodeName,dataAttributes{1,2},dataAttributes{k,1}))
                    end
                    eval(sprintf('vtpRoot(1).%s(1).%s.%s=numArray;',nodeName,dataAttributes{1,2},'data'))
                end
            end

        end
    end

    for i=1:2
        eval(sprintf('vtpRoot(1).Dimensions.%s=str2double(char(pieceDimensions{i}.getValue));',char(pieceDimensions{i}.getName)))
    end
    
    
    
    %vtpRoot.Dimensions.lines=str2double(char(pieceDimensions{1}.getValue));
    %vtpRoot.Dimensions.points=str2double(char(pieceDimensions{2}.getValue));

end
