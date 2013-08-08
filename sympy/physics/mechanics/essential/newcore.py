from essential import ReferenceFrame, Vector, dynamicsymbols
from sympy import diff, integrate, sympify, cacheit, trigsimp, expand, \
     Symbol, solve, sin, acos


def get_motion_pos(position=0, frame=None):

    """
    Calculates the three motion parameters - position, velocity and acceleration
    as vectorial functions of time given the position vector as a function of time.

    Returns a list of vectors as - [acceleration, velocity, position]

    Can also be used for calculations pertaining to rotational motion.

    Parameters
    ==========
    
    position : Vector/VectMul/VectAdd
        Position vector of an object/frame as a function of time

    frame : MovingRefFrame
        The frame to express the motion parameters in

    Examples
    ========

    """
    
    _check_frame(frame)
    if position != 0:
        _check_vector(position)
    else:
        return [0, 0, 0]
    vel = frame.dt(position)
    acc = frame.dt(vel)
    return [acc, vel, position]


def get_motion_vel(velocity=0, position=0, timevalue=0, frame=None):
    """
    Calculates the three motion parameters - position, velocity and acceleration
    as vectorial functions of time given the velocity and a boundary
    condition(position vector at time = timevalue).

    Returns a list of vectors as - [acceleration, velocity, position]

    Can also be used for calculations pertaining to rotational motion.

    Parameters
    ==========

    velocity : Vector/VectMul/VectAdd
        Velocity vector of an object/frame as a function of time

    position : Vector/VectMul/VectAdd
        Boundary condition of position at time = timevalue

    timevalue : Number
        The numeric value of time at which the given boundary condition
        has been expressed

    frame : MovingRefFrame
        The frame to express the motion parameters in

    Examples
    ========

    """

    _check_vector(velocity)
    _check_vector(position)
    timevalue = sympify(timevalue)
    return _process_vector_differential(velocity, position, dynamicsymbols._t,
                                        timevalue, frame)


def get_motion_acc(acceleration=0, velocity=0, position=0, timevalue1=0,
                   timevalue2=0, frame=None):
    """
    Calculates the three motion parameters - position, velocity and acceleration
    as vectorial functions of time given the acceleration and two boundary
    conditions-
    velocity vector at time = timevalue1 and position vector at time = timevalue2.

    Returns a list of [acceleration, velocity, position]

    Parameters
    ==========

    acceleration : Vector/VectMul/VectAdd
        Acceleration of the object/frame as a function of time

    velocity : Vector/VectMul/VectAdd
        Boundary condition of velocity at time = timevalue1

    position : Vector/VectMul/VectAdd
        Boundary condition of position at time = timevalue2

    timevalue1, timevalue2 : Number
        The numeric values of times at which the given boundary conditions
        have been expressed

    frame : MovingRefFrame
        The frame to express the motion parameters in

    Examples
    ========

    """

    _check_vector(acceleration)
    _check_vector(velocity)
    _check_vector(position)
    timevalue1 = sympify(timevalue1)
    timevalue2 = sympify(timevalue2)
    vel = _process_vector_differential(acceleration, velocity, dynamicsymbols._t,
                                       timevalue1, frame)[2]
    pos = _process_vector_differential(vel, position, dynamicsymbols._t,
                                       timevalue2, frame)[2]
    return [acceleration, vel, pos]


def _process_vector_differential(vectdiff, condition, variable, valueofvar, frame):
    """
    Helper function for get_motion methods

    Returns a list of - derivative, vectdiff and integral

    """

    #Make sure boundary condition is independent of 'variable'
    if condition != 0:
        condition = frame.express(condition)
    #Special case of vectdiff == 0
    if vectdiff == 0:
        return [0, 0, condition]
    #Express vectdiff completely in condition's frame to give vectdiff1
    vectdiff1 = frame.express(vectdiff)
    #Find derivative of vectdiff
    vectdiff2 = frame.dt(vectdiff)
    #Integrate and use boundary condition
    vectdiff0 = 0
    for dim in frame:
        function1 = vectdiff1.dot(dim)
        vectdiff0 += _integrate_boundary(function1, variable, valueofvar,
                                         dim.dot(condition)) * dim
    #Return list
    return [vectdiff2, vectdiff, vectdiff0]


def _integrate_boundary(expr, var, valueofvar, value):
    """
    Returns indefinite integratal of expr wrt var, using the boundary
    condition of expr's value being 'value' at var = valueofvar.
    """

    CoI = Symbol('CoI')
    expr = integrate(expr, var) + CoI
    n = expr.subs({CoI : solve(expr.subs({var : valueofvar}) -\
                               value.subs({var : valueofvar}), CoI)[0]})
    return n


def _check_vector(test_vect):
    """
    Helper to check whether an instance is a vector
    """

    if test_vect != 0:
        if not type(test_vect)==Vector:
            raise TypeError(str(test_vect) + " should be a vector.")


class MovingRefFrame(ReferenceFrame):
    def __init__(self, name, pos_vector=None, trans_vel=None, trans_acc=None, \
                 orient_type=None, orient_amount=None, orient_order=None, \
                 ang_vel=None, ang_acc=None, parentframe=None, **kwargs):
        #Special case of no parent frame
        if parentframe is None:
            self._parent = None
            self._root = self
            if 'timevar' not in kwargs:
                #Default time variable
                self._time = dynamicsymbols._t
            else:
                #Set global time
                kwargs['timevar'] = sympify(kwargs['timevar'])
                if not kwargs['timevar'].is_Symbol:
                    raise TypeError("timevar must be a Symbol")
                self._time = kwargs['timevar']
            #All motion params zero
            self._pos_vector = 0
            self._trans_vel = 0
            self._trans_acc = 0
            self._ang_vel = 0
            self._ang_acc = 0
            super(MovingRefFrame, self).__init__(name)
            
        else:
            if not isinstance(parentframe, MovingRefFrame):
                raise TypeError("parentframe should be of type MovingRefFrame")
            self._parent = parentframe
            self._root = parentframe._root
            #Take time variable from parent frame
            self._time = parentframe.time
            #Initial orientation and positioning in case of tuple args
            flag = False
            for x in (trans_vel, trans_acc, ang_vel, ang_acc):
                if type(x) == tuple:
                    flag = True
                    break
            if flag:
                #Time values for boundary conditions must all be zero if tuple args
                #are to be processed. If not, raise exception.
                for x in ('t', 'rt', 't1', 't2', 'rt1', 'rt2'):
                    if x in kwargs:
                        if kwargs[x] != 0:
                            raise ValueError("initial values(t=0) needed for\
                                             initialization with tuple args")
                #Fix initial position from given args/ kwargs
                if pos_vector is not None:
                    kwargs['pos_vector_b'] = pos_vector
                elif 'pos_vector_b' not in kwargs:
                    kwargs['pos_vector_b'] = 0
                self._pos_vector = kwargs['pos_vector_b']
                #Fix initial orientation from given args/ kwargs
                orient_type_temp = 'Axis'
                orient_amount_temp = [0, parentframe.z]
                if orient_type is not None:
                    orient_type_temp = orient_type
                    orient_amount_temp = orient_amount
                elif 'rotation_b' in kwargs:
                    kwargs['rotation_b'] = rotation_vector(tuple(kwargs['rotation_b']), parentframe)
                    orient_type_temp = 'Axis'
                    orient_amount_temp = [kwargs['rotation'].magnitude(),
                                          kwargs['rotation'].normalize()]
                #Calling the superclass' __init__ function for the first time sets
                #the initial position and orientation of this frame wrt parent frame.
                #This allows usage of this frame's basis vectors in setting velocites/
                #acceleration
                super(MovingRefFrame, self).__init__(name)
                super(MovingRefFrame, self).orient(parentframe, orient_type_temp, orient_amount_temp)
                
            #Set translational params as functions of time
            if pos_vector is not None:
                #User has provided pos_vector, hence set motion according to that
                self._trans_acc, self._trans_vel, self._pos_vector = \
                                 get_motion_pos(pos_vector, parentframe)
            elif trans_vel is not None:
                #User has provided trans_vel. Process rest of translational
                #motion params from this
                trans_vel = self._to_vector(trans_vel)
                for x in ('pos_vector_b', 't'):
                    if x not in kwargs:
                        kwargs[x] = 0
                self._trans_acc, self._trans_vel, self._pos_vector = \
                                 get_motion_vel(trans_vel, kwargs['pos_vector_b'],
                                                kwargs['t'], parentframe)
            elif trans_acc is not None:
                #User has provided trans_acc. Process rest of translational motion params
                #using this
                trans_acc = self._to_vector(trans_acc)
                for x in ('trans_vel_b', 'pos_vector_b', 't1', 't2'):
                    if x not in kwargs:
                        kwargs[x] = 0
                self._trans_acc, self._trans_vel, self._pos_vector = \
                                 get_motion_acc(trans_acc, kwargs['trans_vel_b'],
                                                kwargs['pos_vector_b'],
                                                kwargs['t1'], kwargs['t2'], parentframe)
            else:
                #None of the params are provided. This means this frame has no
                #translational motion wrt parent
                self._trans_acc, self._trans_vel, self._pos_vector = (0, 0, 0)
            #Set rotational params as functions of time
            if orient_type is not None or (orient_type is None and ang_vel is None and ang_acc is None):
                super(MovingRefFrame, self).__init__(name)
                if orient_type is None:
                    orient_type = 'Axis'
                    orient_amount = [0, parentframe.z]
                super(MovingRefFrame, self).orient(parentframe, orient_type, orient_amount)
                #Set angular velocity and angular accln params by
                #time-differentiation of DCM
                dcm2diff = self.dcm(parentframe)
                diffed = dcm2diff.diff(self.time)
                angvelmat = diffed * dcm2diff.T
                w1 = trigsimp(expand(angvelmat[7]), recursive=True)
                w2 = trigsimp(expand(angvelmat[2]), recursive=True)
                w3 = trigsimp(expand(angvelmat[3]), recursive=True)
                self._ang_vel = -w1 * parentframe.x - w2 * parentframe.y - \
                                w3 * parentframe.z
                self._ang_acc = parentframe.dt(self._ang_vel)
            elif ang_vel is not None:
                #ang_vel is provided. Process other rotation params using that
                #(and boundary conditions)
                ang_vel = self._to_vector(ang_vel)
                for x in ('rotation_b', 'rt'):
                    if x not in kwargs:
                        kwargs[x] = 0
                self._ang_acc, self._ang_vel, rotation = \
                                 get_motion_vel(ang_vel, \
                                                rotation_vector(tuple(kwargs['rotation_b']), parentframe),
                                                kwargs['rt'], parentframe)
                angle = rotation.magnitude()
                angle.simplify()
                axis = rotation.normalize()
                axis.simplify()
                super(MovingRefFrame, self).__init__(name)
                super(MovingRefFrame, self).orient(parentframe, 'Axis', [angle, axis])
            elif ang_acc is not None:
                #ang_acc is provided. Process other rotation params using that
                #and boundary conditions
                ang_acc = self._to_vector(ang_acc)
                for x in ('ang_vel_b', 'rotation_b', 'rt1', 'rt2'):
                    if x not in kwargs:
                        kwargs[x] = 0
                self._ang_acc, self._ang_vel, rotation = \
                                 get_motion_acc(ang_acc, kwargs['ang_vel_b'],
                                                rotation_vector(tuple(kwargs['rotation_b']), parentframe),
                                                kwargs['rt1'], kwargs['rt2'],
                                                parentframe)
                angle = rotation.magnitude()
                angle.simplify()
                axis = rotation.normalize()
                axis.simplify()
                super(MovingRefFrame, self).__init__(name)
                super(MovingRefFrame, self).orient(parentframe, 'Axis', [angle, axis])
    
    @property
    def parent(self):
        return self._parent

    @property
    def time(self):
        return self._time

    def _to_vector(self, inputparam):
        """
        Converts input given by the user into a vector.

        Input may be a vector, or a tuple defining a vector using two components -
        one in own's basis vectors, and the other in some other frame
        """
        if type(inputparam) != tuple:
            return inputparam
        else:
            #Process tuple to get a vector
            outvect = 0
            outvect += inputparam[0][0] * self.x
            outvect += inputparam[0][1] * self.y
            outvect += inputparam[0][2] * self.z
            outvect += inputparam[1]
            return outvect
        
    def _frame_path(self, otherframe):
        """
        Calculates 'path' of frames starting from this frame to the other,
        along with the index of the common root

        Returns index, list pair
        """
        if self._root != otherframe._root:
            raise ValueError("No connecting path between the two frames-"+ \
                             str(self) + " and " + str(otherframe))
        other_path = []
        frame = otherframe
        while frame.parent is not None:
            other_path.append(frame)
            frame = frame.parent
        other_path.append(frame)
        frameset = set(other_path)
        self_path = []
        frame = self
        while frame not in frameset:
            self_path.append(frame)
            frame = frame.parent
        index = len(self_path)
        i = other_path.index(frame)
        while i >= 0:
            self_path.append(other_path[i])
            i -= 1
        return index, self_path
            
            
    def convert_pos_vector(self, pos_vector, frame):
        """
        Convert a position vector defined in another frame to this frame

        Parameters
        ==========

        pos_vector : vector
            The position vector to be converted

        frame : MovingRefFrame
            The frame the given position vector is defined in

        Examples
        ========

        >>> from sympy.physics.mechanics import MovingRefFrame
        >>> R1 = MovingRefFrame('R1', parentframe=None)
        >>> R2 = MovingRefFrame('R2', parentframe=R1, pos_vector = \
                                2 * R1.basis(0), ang_vel = R1.basis(2))
        >>> pos_vector = 5 * R1.basis(0) + 3 * R1.basis(1) + 4 * R1.basis(2)
        >>> R2.convert_pos_vector(R1)
        ...
        """
        
        
        return pos_vector + frame.pos_vector_in(self)

    
    def pos_vector_in(self, otherframe):
        """
        Returns the relative position vector of this frame's origin in
        another frame.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the position vector in

        Examples
        ========

        >>> from sympy.physics.mechanics import MovingRefFrame
        >>> R1 = MovingRefFrame('R1', parentframe=None)
        >>> R2 = MovingRefFrame('R2', parentframe=R1, ang_vel = R1.basis(2))
        >>> R3 = MovingRefFrame('R3', parentframe = R2, pos_vector = R2.basis(0))
        >>> R3.pos_vector_in(R1)
        #Something equivalent to R2.x
        """

        if otherframe == self:
            return 0
        elif otherframe == self.parent:
            return self._pos_vector
        rootindex, path = self._frame_path(otherframe)
        result = 0
        i = -1
        for i in range(rootindex):
            result += path[i]._pos_vector
        i += 2
        while i < len(path):
            result -= path[i]._pos_vector
            i += 1
        return result
        
    @cacheit
    def trans_vel_in(self, otherframe):
        """
        Returns the relative translational velocity vector of this frame's
        origin in another frame.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative velocity in

        Examples
        ========

        >>> from sympy.physics.mechanics import MovingRefFrame
        >>> R1 = MovingRefFrame('R1', parentframe=None)
        >>> R2 = MovingRefFrame('R2', parentframe=R1, ang_vel = R1.basis(2))
        >>> R3 = MovingRefFrame('R3', parentframe = R2, pos_vector = R2.basis(0))
        >>> R3.trans_vel_in(R1)
        #Something equivalent to R2.y
        """

        if otherframe == self:
            return 0
        elif otherframe == self.parent:
            return self._trans_vel
        elif otherframe.parent == self:
            return -1 * otherframe._trans_vel
        return otherframe.dt(self.pos_vector_in(otherframe))
        
    @cacheit
    def trans_acc_in(self, otherframe):
        """
        Returns the relative translational acceleration vector of this frame's
        origin in another frame.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative acceleration in

        Examples
        ========

        >>> from sympy.physics.mechanics import MovingRefFrame
        >>> R1 = MovingRefFrame('R1', parentframe=None)
        >>> R2 = MovingRefFrame('R2', parentframe=R1, ang_vel = R1.basis(2))
        >>> R3 = MovingRefFrame('R3', parentframe = R2, pos_vector = R2.basis(0))
        >>> R3.trans_vel_in(R1)
        #Something equivalent to -R1.x
        """
        
        if otherframe == self:
            return 0
        elif otherframe == self.parent:
            return self._trans_acc
        elif otherframe.parent == self:
            return -1 * otherframe._trans_acc
        return otherframe.dt(self.trans_vel_in(otherframe))
    
    @cacheit
    def ang_vel_in(self, otherframe):
        """
        Returns the relative angular velocity vector of this frame in
        another frame.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative angular velocity in

        Examples
        ========

        ToBeDone
        """
        
        if otherframe == self:
            return 0
        elif otherframe == self.parent:
            return self._ang_vel
        elif otherframe.parent == self:
            return -1 * otherframe._ang_vel
        rootindex, path = self._frame_path(otherframe)
        result = 0
        i = -1
        for i in range(rootindex):
            result += path[i]._ang_vel
        i += 2
        while i < len(path):
            result -= path[i]._ang_vel
            i += 1
        return result
    
    @cacheit
    def ang_acc_in(self, otherframe):
        """
        Returns the relative angular acceleration vector of this frame in
        another frame.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative angular acceleration in

        Examples
        ========

        ToBeDone
        """
        
        if otherframe == self:
            return 0
        elif otherframe == self.parent:
            return self._ang_acc
        elif otherframe.parent == self:
            return -1 * otherframe._ang_acc
        return otherframe.dt(self.ang_vel_in(otherframe))

    def var_dict(self, otherframe, simplify_mapping=False):
        _check_frame(otherframe)
        vars_matrix = self.dcm(otherframe) * Matrix(otherframe.varlist)
        p_v = self.express(otherframe.pos_vector_in(self))
        mapping = {}
        for i, x in enumerate(self):
            mapping[self.varlist[i]] = simplify(vars_matrix[i] + x.dot(p_v))
        return mapping

    def express(self, vector, variables=False):
        """
        Re-express a vector/scalar in this field
        """
        
        if vector == 0:
            return 0
        if type(vector) == Vector:
            if variables:
                frame_list = [x[-1] for x in vector.args]
                subs_dict = {}
                for frame in frame_list:
                    subs_dict.update(frame.var_dict(self))
                vector = vector.subs(subs_dict).express(self)
                return vector
            else:
                return vector.express(self)
        else:
            frame_list = set([])
            vector = sympify(vector)
            for x in vector.atoms():
                if isinstance(x, sympy.physics.mechanics.essential.CoordinateSym):
                    if x.system not in frame_list and x.system != self:
                        frame_list.add(x.system)
            frame_list = list(frame_list)
            subs_dict = {}
            for frame in frame_list:
                subs_dict.update(frame.var_dict(self))
            return vector.subs(subs_dict)

    def dt(self, vector):
        if vector == 0:
            return sympify(0)
        frame_list = [x[-1] for x in vector.args]
        subs_dict = {}
        for frame in frame_list:
            subs_dict.update(frame.var_dict(self))
        vector = vector.subs(subs_dict)
        return vector.dt(self)

    def __getitem__(self, name):
        if name == 0:
            return self.varlist[0]
        elif name == 1:
            return self.varlist[1]
        elif name == 2:
            return self.varlist[2]
        else:
            raise ValueError("Wrong index")

    def dt(self, expr, order=1):
        """
        Calculate the time derivative of a field function in this frame.

        References
        ==========

        http://en.wikipedia.org/wiki/Rotating_reference_frame#Time_derivatives_in_the_two_frames

        Parameters
        ==========

        expr : vector/scalar Expr
            The field whose time derivative is to be calculated

        order : integer
            The order of the derivative to be calculated

        Examples
        ========

        ToBeDone
        """

        t = dynamicsymbols._t
        if order == 0:
            return expr
        if order%1 != 0 or order < 0:
            raise ValueError("Unsupported value of order entered")
        if isinstance(expr, Vector):
            frame_dict = {}
            diff_dir = {}
            for x in expr.args:
                frame_dict[x[1]] = 0
                diff_dir[x[1]] = 0
                for i, y in enumerate(x[1]):
                    frame_dict[x[1]] += x[0][i] * y
                    diff_dir[x[1]] += diff(x[0][i], t) * y
            #Process each constituent separately, and add to get result
            outvect = 0
            for frame in frame_dict:
                if frame == self:
                    outvect += diff_dir[frame]
                else:
                    outvect += diff_dir[frame] + \
                          frame.ang_vel_in(self).cross(frame_dict[frame])
            return self.dt(outvect, order-1)
        else:
            return diff(self.express(expr), self.time, order)


def _check_frame(test_frame):
    """
    Helper to check whether an instance is a MovingRefFrame
    """

    if not isinstance(test_frame, MovingRefFrame):
        raise ValueError(str(test_frame) + " is not an instance of MovingRefFrame")

def rotation_vector(rotation, parent):
    if isinstance(rotation, tuple):
        if len(rotation) == 2:
            orient_order = None
        else:
            orient_order = rotation[2]
        temp = MovingRefFrame('temp', orient_type = rotation[0],
                              orient_amount = rotation[1],
                              orient_order = orient_order, parentframe = parent)
        dcm = temp.dcm(parent)
        angle = acos(0.5 * (dcm[0] + dcm[4] + dcm[8] - 1))
        e1 = (dcm[7] - dcm[5]) / 2*sin(angle)
        e2 = (dcm[2] - dcm[6]) / 2*sin(angle)
        e3 = (dcm[3] - dcm[1]) / 2*sin(angle)
        axis = (e1 * parent.x + e2 * parent.y + e3 * parent.z).normalize()
        return axis * angle
    else:
        return rotation


class Particle(object):
    """
    A particle with fixed mass and no spatial extension.

    Values need to be supplied at initialization only.
    """

    def __init__(self, name, frame, mass=None, pot_energy=0):
        """
        Initializer for Particle class

        Parameters
        ==========

        name : String
            The name of the new Particle

        frame : MovingRefFrame
            The frame to specify the Particle's motion

        mass : sympifyable
            Mass of the new Particle

        pot_energy : sympifyable
            Potential energy possessed by the Particle
        
        """
        self._mass = sympify(mass)
        if frame is None or not isinstance(frame, MovingRefFrame):
            raise ValueError("Valid frame needs to be entered")
        self._frame = frame
        self._pe = pot_energy

    def __str__(self):
        return name

    @property
    def frame(self):
        """ Returns frame associated with this Particle """
        return self._frame

    @property
    def mass(self):
        """ Returns value of particle's mass """
        return self._mass

    def pos_vector_wrt(self, other):
        """
        The position vector of this particle wrt a specified frame
        or Particle.

        Parameters
        ==========

        other : MovingRefFrame/Particle
            The frame/particle wrt which this particle's position vector is
            to be calculated

        Examples
        ========

        """

        if isinstance(other, MovingRefFrame):
            return self.frame.pos_vector_in(other)
        elif isinstance(other, Particle):
            return self.frame.pos_vector_in(other.frame)
        else:
            raise TypeError("Wrong argument type - " + \
                            str(type(other)))

    def trans_vel_wrt(self, other):
        """
        The translational velocity of this particle wrt a specified frame
        or Particle.

        Parameters
        ==========

        other : MovingRefFrame/Particle
            The frame/particle wrt which this particle's translational vel is
            to be calculated

        Examples
        ========

        """

        if isinstance(other, MovingRefFrame):
            return self.frame.trans_vel_in(other)
        elif isinstance(other, Particle):
            return self.frame.trans_vel_in(other.frame)
        else:
            raise TypeError("Wrong type of argument - " + \
                            str(type(other)))

    def trans_acc_wrt(self, frame):
        """
        The translational acceleration of this particle wrt a specified frame
        or Particle.

        Parameters
        ==========

        other : MovingRefFrame/Particle
            The frame/particle wrt which this particle's translational accln is
            to be calculated

        Examples
        ========

        """

        if isinstance(other, MovingRefFrame):
            return self.frame.trans_acc_in(other)
        elif isinstance(other, Particle):
            return self.frame.trans_acc_in(other.frame)
        else:
            raise TypeError("Wrong type of argument - " + \
                            str(type(other)))

    def linear_momentum(self, frame):
        """
        Linear momentum of the particle in the specified frame.

        The linear momentum L, of a particle P, with respect to frame N is
        given by

        L = m * v

        where m is the mass of the particle, and v is the velocity of the
        particle in the frame N.

        Parameters
        ==========

        frame : MovingRefFrame
            The frame to express the linear momentum of this particle in

        Examples
        ========
        
        """

        if self.mass is None:
            raise ValueError("Mass of the particle has not been defined")
        return self.mass * self.trans_vel_wrt(frame)

    def angular_momentum(self, pos_vector, frame):
        """
        Angular momentum of the particle wrt a certain point in the specified
        frame.

        The angular momentum H, about some point O of a particle, P, is given
        by:

        H = r x m * v

        where r is the position vector from point O to the particle P, m is
        the mass of the particle, and v is the velocity of the particle in
        the inertial frame, N.

        Parameters
        ==========

        pos_vector : vector
            Position vector of the point wrt which the angular momentum is to be
            calculated

        frame : MovingRefFrame
            The frame in which the point's position vector is defined

        Examples
        ========

        """

        if self.mass is None:
            raise ValueError("Mass of the particle has not been defined")
        pos_vector = self.frame.convert_pos_vector(pos_vector, frame)
        return (-1 * pos_vector) ^ (self.mass * self.trans_vel_wrt(frame))

    def kinetic_energy(self, frame):
        """Kinetic energy of the particle in the specified frame

        The kinetic energy, T, of a particle, P, is given by

        'T = 1/2 m v^2'

        where m is the mass of particle P, and v is the velocity of the
        particle in the supplied ReferenceFrame.

        Parameters
        ==========

        frame : MovingRefFrame
            The Particle's velocity is typically defined with respect to
            an inertial frame but any relevant frame in which the velocity is
            known can be supplied.

        Examples
        ========

        """

        if self.mass is None:
            raise ValueError("Mass of the particle has not been defined")
        return (self.mass / sympify(2) * self.trans_vel_wrt(frame) &
                self.trans_vel_wrt(frame))

    @property
    def potential_energy(self):
        """ The potential energy of this particle """
        return self._pe

    def total_energy(self, frame):
        """ Total Energy (PE + KE) of this particle """
        return self.kinetic_energy(frame) + self.potential_energy
