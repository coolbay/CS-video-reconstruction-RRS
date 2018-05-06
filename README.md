
This is the matlab implementation of the paper "Zhao, C., Ma, S., Zhang, J., Xiong, R., & Gao, W. (2017). Video compressive sensing reconstruction via reweighted residual sparsity. IEEE Transactions on Circuits and Systems for Video Technology, 27(6), 1182-1195."


## Usage 
1. Uncompress the rar file `dependencies/mh-bcs-spl-1.0-1` into the folder `Phase1_intra/Utilities/`
2. Run `Phase1_intra/Demo.m` to generate the intial recovery for all frames;
3. Run `Phase2_inter/Demo.m` to apply inter RRS to further enhance the quality of nonkey frames.

## Dependencies
1. [mh-bcs-spl-1.0-1](https://drive.google.com/open?id=15QMtsIdGaZnMMhmG-S3ULUyNpFgT8Q9P)


## Cite this work

Please cite the following paper if you use this code. 
```
@article{zhao2017video,
  title={Video compressive sensing reconstruction via reweighted residual sparsity},
  author={Zhao, Chen and Ma, Siwei and Zhang, Jian and Xiong, Ruiqin and Gao, Wen},
  journal={IEEE Transactions on Circuits and Systems for Video Technology},
  volume={27},
  number={6},
  pages={1182--1195},
  year={2017},
  publisher={IEEE}
}
```
