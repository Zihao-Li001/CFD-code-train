        if (iter % 500 == 0) 
        {
            cout<< iter << " loss = " << loss.item<double>() << endl;
 
            std::ofstream filex("results");
            filex << upred << "\n";
            filex.close();
        }
