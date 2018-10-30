function DataGatherer(X, y, Global)
    save(sprintf('[%s]X.txt', Global.identifier), 'X', '-ascii', '-append');
    save(sprintf('[%s]Y.txt', Global.identifier), 'y', '-ascii', '-append');
end