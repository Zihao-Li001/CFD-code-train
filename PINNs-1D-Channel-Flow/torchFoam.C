#include <torch/torch.h>

using namespace std;

class NN
:
    public torch::nn::Module  //NN 
{
    torch::nn::Sequential net_;

public:

    NN()
    {
        net_ = register_module
        (
            "net", 
            torch::nn::Sequential
            (
                torch::nn::Linear(1,8),
                torch::nn::Tanh(),
                torch::nn::Linear(8,1)
            )
        );
    }

    auto forward(torch::Tensor x)
    {
        return net_->forward(x);
    }
};


int main()
{
    int iter = 1;
    int IterationNum = 1000000;
    double learning_rate = 0.001;

    // Loss function and model
    auto crit = torch::nn::MSELoss();
    auto model = std::make_shared<NN>();
    auto opti = 
        std::make_shared<torch::optim::AdamW>
        (
            model->parameters(), 
            torch::optim::AdamWOptions(learning_rate)
        );

    // Computational domain
    double dy = 0.00125;
    auto init = torch::full({}, 1.0);
    auto mesh = 
        torch::arange(-0.05, 0, dy, torch::requires_grad()).unsqueeze(1);
    int cellN = static_cast<int>(0.05/dy);
    auto dpdxByNu = torch::ones({cellN, 1}, torch::kFloat);
    auto nut = torch::zeros({cellN});
        
    for (int i = 0; i < IterationNum; i++) 
    {
        opti->zero_grad();

        // Forward pass
        auto upred = model->forward(mesh);

        // Compute gradients
        auto dudy = 
            torch::autograd::grad
            (
                {upred},
                {mesh},
                {torch::ones_like(upred)},
                true,
                true
            )[0];
        auto dudyy = 
            torch::autograd::grad
            (
                {dudy},
                {mesh},
                {torch::ones_like(upred)},
                true,
                true
            )[0];

        // Boundary conditions and initial conditions
        auto dudyBottom = dudy[0][0];
        auto dudyTop = dudy[cellN - 1][0];
        auto meanU = torch::mean(upred);
        auto ubottom = torch::full_like(upred[0], 0.0375);
        auto dpdxByNu = torch::full({cellN}, 1.2e-05/1e-8);
        //#include "nut40.H"
        //auto dpdxByNu = 0.018/nut;
        //auto ubottom = torch::full_like(upred[0], 0.77);
        
        // Loss computation
        auto loss_bottom = crit(upred[0], ubottom);
        auto loss_top = crit(dudyTop, 0.0*dudyTop);
        auto loss_bd = loss_bottom + 1000*loss_top;
        auto loss_ini = 1000*crit(meanU, init);
        auto loss_pde = crit(-dudyy.reshape({cellN}), dpdxByNu);
        auto loss = loss_pde + loss_bd + loss_ini;
        //auto loss = loss_pde + loss_bd;

        // Backward pass and optimization step
        loss.backward();
        opti->step();
        
        #include "output.H"
        iter++;
    }

    // Save the model
    torch::save(model, "model.pth");
    std::cout<< "Done!" << std::endl;
    return 0;
}
