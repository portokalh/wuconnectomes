


class registrationparams(object):

    def __init__(self, nbins=32, sampling_proportion=None, level_iters=[10000, 1000, 100], ):
        r"""Initialize an instance of the Mutual Information metric.

        This class implements the methods required by Optimizer to drive the
        registration process.

        Parameters
        ----------
        nbins : int, optional
            the number of bins to be used for computing the intensity
            histograms. The default is 32.
        sampling_proportion : None or float in interval (0, 1], optional
            There are two types of sampling: dense and sparse. Dense sampling
            uses all voxels for estimating the (joint and marginal) intensity
            histograms, while sparse sampling uses a subset of them. If
            `sampling_proportion` is None, then dense sampling is
            used. If `sampling_proportion` is a floating point value in (0,1]
            then sparse sampling is used, where `sampling_proportion`
            specifies the proportion of voxels to be used. The default is
            None.

        Notes
        -----
        Since we use linear interpolation, images are not, in general,
        differentiable at exact voxel coordinates, but they are differentiable
        between voxel coordinates. When using sparse sampling, selected voxels
        are slightly moved by adding a small random displacement within one
        voxel to prevent sampling points from being located exactly at voxel
        coordinates. When using dense sampling, this random displacement is
        not applied.

        """
        self.histogram = ParzenJointHistogram(nbins)
        self.sampling_proportion = sampling_proportion
        self.metric_val = None
        self.metric_grad = None