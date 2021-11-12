classdef Eqf_Filter_Robot  < handle
    properties
        %Initial Conditions of the Fiter
        Q_0;
        State_Observer;
        Output_Matrix_Gain;
        Rica_Coef;
        Rica;
        State_Matrix_Gain_Coef;
        State_Matrix_Gain;
        Initial_Estimated_Robot_Pose;
        Robot_State_Observer;
        B_t_Matrix
        %Temperary Data Collected within Filter
        dt;
        Robot_R_t;
        Landmark_R_t;
        Landmark_Quant;
        Est_Landmark_Position;
        De_m_M;
        Output_Matrix_M;
        State_Matrix_M;
        delta_M;
        Landmark_Est;
        Est_Robot_Pose; 
        Body_Fixed_Landmark;
        R_t;
        Q_i;
        Uncorrected_Pose;
    end
    
    methods
        
        function obj = Eqf_Filter_Robot(Landmark_Quant)
            obj.Landmark_Quant = Landmark_Quant;
        end
        
        function Time_Change(obj,dt)
            obj.dt = dt;
            if dt < 0
                 warning("Negative Time Step");
            end
        end
        
        function Filter_Initiation(obj,Landmark_Quant)
            obj.Landmark_Quant = Landmark_Quant;
            %without velocity bias
            obj.Rica = obj.Rica_Coef * eye(3 * Landmark_Quant + 6);
            obj.State_Matrix_Gain = obj.State_Matrix_Gain_Coef * eye(3 * Landmark_Quant + 6);
            %with velocity bias
%             obj.Rica = obj.Rica_Coef * eye(3 * Landmark_Quant + 6);
%             obj.State_Matrix_Gain = obj.State_Matrix_Gain_Coef * eye(3 * Landmark_Quant + 6);
        end
        
        function Inertial_Initiate(obj,Landmarks_Initial_Est, Est_Robot_Pose)
            Q_0 = [];
            State_Observer = [];
            for i = 0:obj.Landmark_Quant - 1
                State_Observer = [State_Observer eye(4)];
                Q_0_i = Est_Robot_Pose(1:3,1:3)' * (Landmarks_Initial_Est(1,i * 3 + 1: (i + 1) * 3)' ...
                    - Est_Robot_Pose(1:3,4));
                Q_0 = [Q_0 Q_0_i'];
            end
            obj.Q_0 = Q_0;
            obj.State_Observer = State_Observer;
        end
        
         function Body_Initiate(obj,Landmarks_Initial_Est)
            Q_0 = [];
            State_Observer = [];
            for i = 0:obj.Landmark_Quant - 1
                State_Observer = [State_Observer eye(4)];
                Q_0_i = Landmarks_Initial_Est(1,i * 3 + 1: (i + 1) * 3)';
                Q_0 = [Q_0 Q_0_i'];
            end
            obj.Q_0 = Q_0;
            obj.State_Observer = State_Observer;
         end
         
         function Robot_Pose_Initiate(obj,Est_Robot_Pose)
             obj.Initial_Estimated_Robot_Pose = Est_Robot_Pose;
             obj.Robot_State_Observer = eye(4);
             obj.Uncorrected_Pose = Est_Robot_Pose;
         end
         
         function  Motion_Update(obj, Omega, Velocity, Frame)
             Landmark_Est = [];
             Body_Fixed_Landmark_M = [];
             obj.State_Matrix_M = eye(1 * obj.Landmark_Quant);
             obj.Output_Matrix_M = zeros(obj.Landmark_Quant, 3*obj.Landmark_Quant);
             B_t = []; 
             Q_i_M = [];
             %The robot state gets updated by the first term of the
             %innovation here 
             U = [obj.skew(Omega) Velocity;
                  zeros(1,4)];
             %Do we use Est_Robot_Pose or State_Observer here 
             obj.Robot_State_Observer = obj.Robot_State_Observer * expm(obj.dt * U);
             obj.Est_Robot_Pose = obj.Initial_Estimated_Robot_Pose * obj.Robot_State_Observer;
             obj.Uncorrected_Pose = obj.Uncorrected_Pose * expm(obj.dt * U);
             %Do we use Est_Robot_Pose or State_Observer here 
             Robot_B_t = [obj.Robot_State_Observer(1:3,1:3),zeros(3,3);
                          obj.skew(obj.Robot_State_Observer(1:3,4)) * obj.Robot_State_Observer(1:3,1:3) , obj.Robot_State_Observer(1:3,1:3)];
             for j = 0:obj.Landmark_Quant - 1
                %Basic Motion Estimation
                State_Observer_i = obj.State_Observer(1:4, j * 4 + 1:(j+1) * 4);
                State_Observer_R = State_Observer_i(1:3,1:3);
                State_Observer_C = State_Observer_i(4,4);
                
                %Single Landmark Estimation
                Initial_Position = obj.Q_0(1,j * 3 + 1:(j+1) * 3)';
                Estimated_Position = (1/State_Observer_C) * State_Observer_R' * Initial_Position;
                Body_Fixed_Landmark_M = [Body_Fixed_Landmark_M Estimated_Position'];
                
                %Velocity Noise
                Lambda_c = Estimated_Position' * Velocity / (Estimated_Position'*Estimated_Position);
                Lambda_R = Omega + cross(Estimated_Position,Velocity)/ (Estimated_Position' * Estimated_Position);
                Lambda = [obj.skew(Lambda_R) zeros(3,1);
                          zeros(1,3)         Lambda_c];
                State_Observer_i = State_Observer_i * expm(obj.dt * Lambda);
                
                %V0
                Origin_Velocity  = State_Observer_C * State_Observer_R * Velocity;
                
                %State matrix A_t^0 design
                State_Matrix = (Origin_Velocity * Initial_Position' - obj.skew(Origin_Velocity) * obj.skew(Initial_Position))/(Initial_Position' * Initial_Position);
                obj.State_Matrix_M(j*3+1:(j+1)*3,j*3+1:(j+1)*3) = State_Matrix;
                
                %Output matrix C^0 design
                Output_Matrix = (Initial_Position/norm(Initial_Position))';
                obj.Output_Matrix_M(j+1,j*3+1:(j+1)*3) = Output_Matrix;
                
                %Velocity Bias B_t design
                Q_i = obj.Est_Robot_Pose(1:3,1:3)' * (Estimated_Position - obj.Est_Robot_Pose(1:3,4));
                Q_i_M = [Q_i_M Q_i]; 
                %Do we use Est_Robot_Pose or State_Observer here 
                B_t_i_V = obj.skew(Initial_Position) * obj.State_Observer(1:3,1:3) * (obj.skew(Q_i)/ (Q_i' * Q_i)) ...
                    - Initial_Position * Q_i' / (Q_i' * Q_i);
                B_t_i_Omega = obj.skew(Initial_Position); 
                B_t_i = [B_t_i_Omega B_t_i_V];
                %B_t_i = [B_t_i_V B_t_i_Omega];
                B_t = [B_t; B_t_i];
                %Transfer Body Frame to Inertial Frame
                if strcmp(Frame, 'Inertial')
                    Estimated_Position = obj.Est_Robot_Pose(1:3,1:3) * Estimated_Position + obj.Est_Robot_Pose(1:3,4);
                end
                %This needs initial pose
                obj.Body_Fixed_Landmark = Body_Fixed_Landmark_M;
                Landmark_Est = [Landmark_Est Estimated_Position'];
                obj.Landmark_Est = Landmark_Est;
                obj.State_Observer(1:4, j * 4 + 1:(j+1) * 4) = State_Observer_i;
             end
            obj.Q_i = Q_i_M;
            %Ricatti Equation Update
            if min(eig(obj.Rica)) < 0
                 warning("Ricatti Matrix getting too small");
            end
             %Without Velocity Bias
             A_discrete = obj.State_Matrix_M * obj.dt + eye(size(obj.State_Matrix_M));
%             obj.Rica = A_discrete * obj.Rica * A_discrete' + obj.dt * obj.State_Matrix_Gain;
%             if min(eig(obj.Rica)) < 0
%                  warning("Ricatti Matrix getting too small");
%             end
            obj.State_Matrix_M = [zeros(6,18); zeros(12,6) obj.State_Matrix_M];
            obj.Output_Matrix_M = [zeros(4,6) obj.Output_Matrix_M];
            B_t_Matrix = [Robot_B_t; B_t];
            obj.B_t_Matrix = B_t_Matrix * obj.Robot_R_t * B_t_Matrix';
            %obj.State_Matrix_Gain(1:6,1:6) = 0;
            Rica_dot = obj.State_Matrix_M * obj.Rica + obj.Rica * obj.State_Matrix_M'  + obj.B_t_Matrix + obj.State_Matrix_Gain;
            obj.Rica = obj.Rica + Rica_dot * obj.dt;
            obj.delta_M = zeros(obj.Landmark_Quant,1);
            obj.De_m_M = zeros(4*obj.Landmark_Quant,3*obj.Landmark_Quant);
        end
        
        function Delta_Calculation(obj,Range_Measurement,Index)
            j = Index;
            %Calling the variables
            State_Observer_i = obj.State_Observer(1:4, j * 4 + 1:(j+1) * 4);
            State_Observer_C = State_Observer_i(4,4);
            
            % Create the measurement error vector delta
            Error_Measurement = State_Observer_C * Range_Measurement;
            Initial_Position = obj.Q_0(1,j * 3 + 1:(j+1) * 3)';
            delta = Error_Measurement - norm(Initial_Position);
            %A sample of noise rejection module
              if delta == -norm(Initial_Position)  %|| delta < -0.3
                  delta = 0;
              end
            obj.delta_M(j+1) = delta;
            
            %DE|id (first term of innovation)
            De_m = [-obj.skew(Initial_Position)/(Initial_Position'*Initial_Position);
                -Initial_Position'/(Initial_Position'*Initial_Position)];
            obj.De_m_M(j*4+1:(j+1)*4,j*3+1:(j+1)*3) = De_m;
        end
        function Innovation(obj, Landmarks_Initial_True)
            Innovation_Update = obj.Rica * obj.Output_Matrix_M' * obj.Output_Matrix_Gain^(-1) * obj.delta_M;
            Innovation = obj.De_m_M * Innovation_Update(7:end,1);
            Robot_Innovation = ([obj.skew(Innovation_Update(1:3)), Innovation_Update(4:6);
                                    zeros(1,4)]);
            obj.Robot_State_Observer = expm(obj.dt * Robot_Innovation) * obj.Robot_State_Observer;
            %disp(obj.Initial_Estimated_Robot_Pose * obj.Robot_State_Observer);
            %disp(obj.Uncorrected_Pose);
            %Innovation = obj.De_m_M * obj.Rica * [zeros(4,6) obj.Output_Matrix_M]' * obj.Output_Matrix_Gain^(-1) * obj.delta_M;
            %Innovation = obj.De_m_M * obj.Rica(7:end,7:end) * obj.Output_Matrix_M' * obj.Output_Matrix_Gain^(-1) * obj.delta_M;
            Eqf_M = [];
            for j = 0:obj.Landmark_Quant-1
                %innovation update
                Delta_R = Innovation(4*j + (1:3));
                Delta_C = Innovation(4*j + 4);
                Delta = [obj.skew(Delta_R) zeros(3,1);
                         zeros(1,3)        Delta_C];
                
                %Basic Motion Estimation
                State_Observer_i = obj.State_Observer(1:4, j * 4 + 1:(j+1) * 4);
                State_Observer_R = State_Observer_i(1:3,1:3);
                State_Observer_C = State_Observer_i(4,4);
                
                %State and equivariant filter errror estimation
                State_Observer_i = expm(obj.dt * Delta) * State_Observer_i;
                %Change to body fixed frame
                eqf_e = State_Observer_C * State_Observer_R * Landmarks_Initial_True(1,j * 3 + 1:(j+1) * 3)';
                %eqf_e = (Est_Rotation')^-1 * eqf_e + Assumed_Position;
                %Eqf_Error = [Eqf_Error eqf_e'];
                obj.State_Observer(1:4, j * 4 + 1:(j+1) * 4) = State_Observer_i;
                Eqf_M = [Eqf_M eqf_e'];
            end
            if min(eig(obj.Rica)) < 0
                 warning("Ricatti Matrix getting too small");
            end
            %without velocity bias
           % obj.Rica = (obj.dt * obj.Output_Matrix_M' * obj.Output_Matrix_Gain^-1 * obj.Output_Matrix_M + obj.Rica^(-1))^(-1);
            %             Rica_dot = - obj.Rica * obj.Output_Matrix_M' * obj.Output_Matrix_Gain^-1 * obj.Output_Matrix_M * obj.Rica;
            %with velocity bias
             Rica_dot = - obj.Rica *  obj.Output_Matrix_M' * obj.Output_Matrix_Gain^-1 * obj.Output_Matrix_M * obj.Rica;
             obj.Rica = obj.Rica + Rica_dot * obj.dt;
            if min(eig(obj.Rica)) < 0
                warning("Ricatti Matrix getting too small");
             end
        end
        function m = skew(obj,v)
            %creating check functions 
            assert(numel(v) == 3, "v is not a 3-vector");
            %creating 3 dimensional skewed matrix
            m = [0   -v(3)  v(2);
                 v(3) 0    -v(1);
                -v(2) v(1)  0];
        end
        
    end
end