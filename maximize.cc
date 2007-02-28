#include <iostream>
#include <bernstein/bernstein.h>
#include <barvinok/evalue.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include <barvinok/bernstein.h>
#include "argp.h"
#include "evalue_convert.h"

using std::cout;
using std::cerr;
using std::endl;
using namespace GiNaC;
using namespace bernstein;
using namespace barvinok;

#define OPT_VARS  	    (BV_OPT_LAST+1)
#define OPT_SPLIT  	    (BV_OPT_LAST+2)
#define OPT_MIN  	    (BV_OPT_LAST+3)

struct argp_option argp_options[] = {
    { "split",		    OPT_SPLIT,	"int" },
    { "variables",	    OPT_VARS,  	"list",	0,
	"comma separated list of variables over which to maximize" },
    { "verbose",	    'V',  	0,	0, },
    { "minimize",	    OPT_MIN,  	0, 0,	"minimize instead of maximize"},
    { 0 }
};

struct options {
    struct convert_options   convert;
    char* var_list;
    int verbose;
    int split;
    int minimize;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct options *options = (struct options*) state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = &options->convert;
	options->var_list = NULL;
	options->verbose = 0;
	options->split = 0;
	options->minimize = 0;
	break;
    case 'V':
	options->verbose = 1;
	break;
    case OPT_VARS:
	options->var_list = strdup(arg);
	break;
    case OPT_SPLIT:
	options->split = atoi(arg);
	break;
    case OPT_MIN:
	options->minimize = 1;
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}


#define ALLOC(type) (type*)malloc(sizeof(type))
#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

enum token_type { TOKEN_UNKNOWN = 256, TOKEN_VALUE, TOKEN_IDENT, TOKEN_GE,
		  TOKEN_UNION };

struct token {
    enum token_type  type;

    int line;
    int col;

    union {
	Value	v;
	char	*s;
    } u;
};

static struct token *token_new(int line, int col)
{
    struct token *tok = ALLOC(struct token);
    tok->line = line;
    tok->col = col;
    return tok;
}

void token_free(struct token *tok)
{
    if (tok->type == TOKEN_VALUE)
	value_clear(tok->u.v);
    else if (tok->type == TOKEN_IDENT)
	free(tok->u.s);
    free(tok);
}

struct stream {
    FILE    	    *file;
    int		    line;
    int		    col;
    int		    eof;

    char	    *buffer;
    size_t	    size;
    size_t	    len;
    int		    c;

    struct token    *tokens[5];
    int		    n_token;
};

static struct stream* stream_new(FILE *file)
{
    int i;
    struct stream *s = ALLOC(struct stream);
    s->file = file;
    s->size = 256;
    s->buffer = (char*)malloc(s->size);
    s->len = 0;
    s->line = 1;
    s->col = 0;
    s->eof = 0;
    s->c = -1;
    for (i = 0; i < 5; ++i)
	s->tokens[i] = NULL;
    s->n_token = 0;
    return s;
}

static void stream_free(struct stream *s)
{
    free(s->buffer);
    assert(s->n_token == 0);
    free(s);
}

static int stream_getc(struct stream *s)
{
    int c;
    if (s->eof)
	return -1;
    c = fgetc(s->file);
    if (c == -1)
	s->eof = 1;
    if (s->c != -1) {
	if (s->c == '\n') {
	    s->line++;
	    s->col = 0;
	} else
	    s->col++;
    }
    s->c = c;
    return c;
}

static void stream_ungetc(struct stream *s, int c)
{
    ungetc(c, s->file);
    s->c = -1;
}

static void stream_push_char(struct stream *s, int c)
{
    if (s->len >= s->size) {
	s->size = (3*s->size)/2;
	s->buffer = (char*)realloc(s->buffer, s->size);
    }
    s->buffer[s->len++] = c;
}

static struct token *stream_push_token(struct stream *s, struct token *tok)
{
    assert(s->n_token < 5);
    s->tokens[s->n_token++] = tok;
}

static struct token *stream_next_token(struct stream *s)
{
    int c;
    struct token *tok;
    int line, col;

    if (s->n_token)
	return s->tokens[--s->n_token];

    s->len = 0;

    /* skip spaces */
    while ((c = stream_getc(s)) != -1 && isspace(c))
	/* nothing */
	;

    line = s->line;
    col = s->col;

    if (c == -1)
	return NULL;
    if (c == '(' ||
        c == ')' ||
	c == '+' ||
	c == '/' ||
	c == '*' ||
	c == '^' ||
	c == '=' ||
	c == ',' ||
	c == '_' ||
	c == '[' ||
	c == ']' ||
	c == '{' ||
	c == '}') {
	tok = token_new(line, col);
	tok->type = (token_type)c;
	return tok;
    }
    if (c == '-' || isdigit(c)) {
	tok = token_new(line, col);
	tok->type = TOKEN_VALUE;
	value_init(tok->u.v);
	stream_push_char(s, c);
	while ((c = stream_getc(s)) != -1 && isdigit(c))
	    stream_push_char(s, c);
	if (c != -1)
	    stream_ungetc(s, c);
	if (s->len == 1 && s->buffer[0] == '-')
	    value_set_si(tok->u.v, -1);
	else {
	    stream_push_char(s, '\0');
	    mpz_set_str(tok->u.v, s->buffer, 0);
	}
	return tok;
    }
    if (isalpha(c)) {
	tok = token_new(line, col);
	stream_push_char(s, c);
	while ((c = stream_getc(s)) != -1 && isalpha(c))
	    stream_push_char(s, c);
	if (c != -1)
	    stream_ungetc(s, c);
	stream_push_char(s, '\0');
	if (!strcmp(s->buffer, "UNION")) {
	    tok->type = TOKEN_UNION;
	} else {
	    tok->type = TOKEN_IDENT;
	    tok->u.s = strdup(s->buffer);
	}
	return tok;
    }
    if (c == '>') {
	if ((c = stream_getc(s)) == '=') {
	    tok = token_new(line, col);
	    tok->type = TOKEN_GE;
	    return tok;
	}
	if (c != -1)
	    stream_ungetc(s, c);
    }

    tok = token_new(line, col);
    tok->type = TOKEN_UNKNOWN;
    return tok;
}

void stream_error(struct stream *s, struct token *tok, char *msg)
{
    int line = tok ? tok->line : s->line;
    int col = tok ? tok->col : s->col;
    fprintf(stderr, "syntax error (%d, %d): %s\n", line, col, msg);
    if (tok) {
	if (tok->type < 256)
	    fprintf(stderr, "got '%c'\n", tok->type);
    }
}

struct parameter {
    char    	    	*name;
    int	     		 pos;
    struct parameter	*next;
};

struct parameter *parameter_new(char *name, int pos, struct parameter *next)
{
    struct parameter *p = ALLOC(struct parameter);
    p->name = strdup(name);
    p->pos = pos;
    p->next = next;
    return p;
}

static int parameter_pos(struct parameter **p, char *s)
{
    int pos = *p ? (*p)->pos+1 : 0;
    struct parameter *q;

    for (q = *p; q; q = q->next) {
	if (strcmp(q->name, s) == 0)
	    break;
    }
    if (q)
	pos = q->pos;
    else
	*p = parameter_new(s, pos, *p);
    return pos;
}

static int optional_power(struct stream *s)
{
    int pow;
    struct token *tok;
    tok = stream_next_token(s);
    if (!tok)
	return 1;
    if (tok->type != '^') {
	stream_push_token(s, tok);
	return 1;
    }
    token_free(tok);
    tok = stream_next_token(s);
    if (!tok || tok->type != TOKEN_VALUE) {
	stream_error(s, tok, "expecting exponent");
	if (tok)
	    stream_push_token(s, tok);
	return 1;
    }
    pow = VALUE_TO_INT(tok->u.v);
    token_free(tok);
    return pow;
}

static evalue *evalue_read_factor(struct stream *s, struct parameter **p);
static evalue *evalue_read_term(struct stream *s, struct parameter **p);

static evalue *create_fract_like(struct stream *s, evalue *arg, enode_type type,
			         struct parameter **p)
{
    evalue *e;
    int pow;
    pow = optional_power(s);

    e = ALLOC(evalue);
    value_init(e->d);
    e->x.p = new_enode(type, pow+2, -1);
    value_clear(e->x.p->arr[0].d);
    e->x.p->arr[0] = *arg;
    free(arg);
    evalue_set_si(&e->x.p->arr[1+pow], 1, 1);
    while (--pow >= 0)
	evalue_set_si(&e->x.p->arr[1+pow], 0, 1);

    return e;
}

static evalue *read_fract(struct stream *s, struct token *tok, struct parameter **p)
{
    evalue *arg;

    tok = stream_next_token(s);
    assert(tok);
    assert(tok->type == '{');

    token_free(tok);
    arg = evalue_read_term(s, p);
    tok = stream_next_token(s);
    if (!tok || tok->type != '}') {
	stream_error(s, tok, "expecting \"}\"");
	if (tok)
	    stream_push_token(s, tok);
    } else
	token_free(tok);

    return create_fract_like(s, arg, fractional, p);
}

static evalue *read_periodic(struct stream *s, struct parameter **p)
{
    evalue **list;
    int len;
    int n;
    evalue *e = NULL;

    struct token *tok;
    tok = stream_next_token(s);
    assert(tok && tok->type == '[');
    token_free(tok);

    len = 100;
    list = (evalue **)malloc(len * sizeof(evalue *));
    n = 0;

    for (;;) {
	evalue *e = evalue_read_term(s, p);
	if (!e) {
	    stream_error(s, NULL, "missing argument or list element");
	    goto out;
	}
	if (n >= len) {
	    len = (3*len)/2;
	    list = (evalue **)realloc(list, len * sizeof(evalue *));
	}
	list[n++] = e;

	tok = stream_next_token(s);
	if (!tok) {
	    stream_error(s, NULL, "unexpected EOF");
	    goto out;
	}
	if (tok->type != ',')
	    break;
	token_free(tok);
    }

    if (tok->type != ']') {
	stream_error(s, tok, "expecting \"]\"");
	stream_push_token(s, tok);
	goto out;
    }

    token_free(tok);

    tok = stream_next_token(s);
    if (tok->type == '_') {
	int pos;
	token_free(tok);
	tok = stream_next_token(s);
	if (!tok || tok->type != TOKEN_IDENT) {
	    stream_error(s, tok, "expecting identifier");
	    if (tok)
		stream_push_token(s, tok);
	    goto out;
	}
	e = ALLOC(evalue);
	value_init(e->d);
	pos = parameter_pos(p, tok->u.s);
	token_free(tok);
	e->x.p = new_enode(periodic, n, pos+1);
	while (--n >= 0) {
	    value_clear(e->x.p->arr[n].d);
	    e->x.p->arr[n] = *list[n];
	    free(list[n]);
	}
    } else if (n == 1) {
	stream_push_token(s, tok);
	e = create_fract_like(s, list[0], flooring, p);
	n = 0;
    } else {
	stream_error(s, tok, "unexpected token");
	stream_push_token(s, tok);
    }

out:
    while (--n >= 0) {
	free_evalue_refs(list[n]);
	free(list[n]);
    }
    free(list);
    return e;
}

static evalue *evalue_read_factor(struct stream *s, struct parameter **p)
{
    struct token *tok;
    evalue *e = NULL;

    tok = stream_next_token(s);
    if (!tok)
	return NULL;

    if (tok->type == '(') {
	token_free(tok);
	e = evalue_read_term(s, p);
	tok = stream_next_token(s);
	if (!tok || tok->type != ')') {
	    stream_error(s, tok, "expecting \")\"");
	    if (tok)
		stream_push_token(s, tok);
	} else
	    token_free(tok);
    } else if (tok->type == TOKEN_VALUE) {
	e = ALLOC(evalue);
	value_init(e->d);
	value_set_si(e->d, 1);
	value_init(e->x.n);
	value_assign(e->x.n, tok->u.v);
	token_free(tok);
	tok = stream_next_token(s);
	if (tok && tok->type == '/') {
	    token_free(tok);
	    tok = stream_next_token(s);
	    if (!tok || tok->type != TOKEN_VALUE) {
		stream_error(s, tok, "expecting denominator");
		if (tok)
		    stream_push_token(s, tok);
		return NULL;
	    }
	    value_assign(e->d, tok->u.v);
	    token_free(tok);
	} else if (tok)
	    stream_push_token(s, tok);
    } else if (tok->type == TOKEN_IDENT) {
	int pos = parameter_pos(p, tok->u.s);
	int pow = optional_power(s);
	token_free(tok);
	e = ALLOC(evalue);
	value_init(e->d);
	e->x.p = new_enode(polynomial, pow+1, pos+1);
	evalue_set_si(&e->x.p->arr[pow], 1, 1);
	while (--pow >= 0)
	    evalue_set_si(&e->x.p->arr[pow], 0, 1);
    } else if (tok->type == '[') {
	stream_push_token(s, tok);
	e = read_periodic(s, p);
    } else if (tok->type == '{') {
	stream_push_token(s, tok);
	e = read_fract(s, tok, p);
    }

    tok = stream_next_token(s);
    if (tok && tok->type == '*') {
	evalue *e2;
	token_free(tok);
	e2 = evalue_read_factor(s, p);
	if (!e2) {
	    stream_error(s, NULL, "unexpected EOF");
	    return NULL;
	}
	emul(e2, e);
	free_evalue_refs(e2);
	free(e2);
    } else if (tok)
	stream_push_token(s, tok);

    return e;
}

static evalue *evalue_read_term(struct stream *s, struct parameter **p)
{
    struct token *tok;
    evalue *e = NULL;

    e = evalue_read_factor(s, p);
    if (!e)
	return NULL;

    tok = stream_next_token(s);
    if (!tok)
	return e;

    if (tok->type == '+') {
	evalue *e2;
	token_free(tok);
	e2 = evalue_read_term(s, p);
	if (!e2) {
	    stream_error(s, NULL, "unexpected EOF");
	    return NULL;
	}
	eadd(e2, e);
	free_evalue_refs(e2);
	free(e2);
    } else
	stream_push_token(s, tok);

    return e;
}

struct constraint {
    int	    		 type;
    Vector  		*v;
    struct constraint 	*next;
    struct constraint 	*union_next;
};

static struct constraint *constraint_new()
{
    struct constraint *c = ALLOC(struct constraint);
    c->type = -1;
    c->v = Vector_Alloc(16);
    c->next = NULL;
    c->union_next = NULL;
    return c;
}

static void constraint_free(struct constraint *c)
{
    while (c) {
	struct constraint *next = c->next ? c->next : c->union_next;
	Vector_Free(c->v);
	free(c);
	c = next;
    }
}

static void constraint_extend(struct constraint *c, int pos)
{
    Vector *v;
    if (pos < c->v->Size)
	return;

    v = Vector_Alloc((3*c->v->Size)/2);
    Vector_Copy(c->v->p, v->p, c->v->Size);
    Vector_Free(c->v);
    c->v = v;
}

static int evalue_read_constraint(struct stream *s, struct parameter **p,
				  struct constraint **constraints,
				  struct constraint *union_next)
{
    struct token *tok;
    struct constraint *c = NULL;

    while ((tok = stream_next_token(s))) {
	struct token *tok2;
	int pos;
	if (tok->type == '+')
	    token_free(tok);
	else if (tok->type == TOKEN_IDENT) {
	    if (!c)
		c = constraint_new();
	    pos = parameter_pos(p, tok->u.s);
	    constraint_extend(c, 1+pos);
	    value_set_si(c->v->p[1+pos], 1);
	    token_free(tok);
	} else if (tok->type == TOKEN_VALUE) {
	    if (!c)
		c = constraint_new();
	    tok2 = stream_next_token(s);
	    if (tok2 && tok2->type == TOKEN_IDENT) {
		pos = parameter_pos(p, tok2->u.s);
		constraint_extend(c, 1+pos);
		value_assign(c->v->p[1+pos], tok->u.v);
		token_free(tok);
		token_free(tok2);
	    } else {
		if (tok2)
		    stream_push_token(s, tok2);
		value_assign(c->v->p[0], tok->u.v);
		token_free(tok);
	    }
	} else if (tok->type == TOKEN_GE || tok->type == '=') {
	    int type = tok->type == TOKEN_GE;
	    token_free(tok);
	    tok = stream_next_token(s);
	    if (!tok || tok->type != TOKEN_VALUE || value_notzero_p(tok->u.v)) {
		stream_error(s, tok, "expecting \"0\"");
		if (tok)
		    stream_push_token(s, tok);
		*constraints = NULL;
	    } else {
		c->type = type;
		c->next = *constraints;
		c->union_next = union_next;
		*constraints = c;
		token_free(tok);
	    }
	    break;
	} else {
	    if (!c)
		stream_push_token(s, tok);
	    else {
		stream_error(s, tok, "unexpected token");
		*constraints = NULL;
	    }
	    return 0;
	}
    }
    return tok != NULL;
}

static struct constraint *evalue_read_domain(struct stream *s, struct parameter **p,
					     unsigned MaxRays)
{
    struct constraint *constraints = NULL;
    struct constraint *union_next = NULL;
    struct token *tok;
    int line;

    tok = stream_next_token(s);
    if (!tok)
	return NULL;
    stream_push_token(s, tok);

    line = tok->line;
    while (evalue_read_constraint(s, p, &constraints, union_next)) {
	tok = stream_next_token(s);
	if (tok) {
	    if (tok->type == TOKEN_UNION) {
		token_free(tok);
		tok = stream_next_token(s);
		if (!tok) {
		    stream_error(s, NULL, "unexpected EOF");
		    return constraints;
		}
		stream_push_token(s, tok);
		union_next = constraints;
		constraints = NULL;
	    } else {
		union_next = NULL;
		stream_push_token(s, tok);
		/* empty line separates domain from evalue */
		if (tok->line > line+1)
		    break;
	    }
	    line = tok->line;
	}
    }
    return constraints;
}

struct section {
    struct constraint	*constraints;
    evalue		*e;
    struct section	*next;
};

static char **extract_parameters(struct parameter *p, unsigned *nparam)
{
    int i;
    char **params;

    *nparam = p ? p->pos+1 : 0;
    params = ALLOCN(char *, *nparam);
    for (i = 0; i < *nparam; ++i) {
	struct parameter *next = p->next;
	params[p->pos] = p->name;
	free(p);
	p = next;
    }
    return params;
}

static Polyhedron *constraints2domain(struct constraint *constraints,
				      unsigned nparam, unsigned MaxRays)
{
    Polyhedron *D;
    Matrix *M;
    int n;
    struct constraint *c;
    struct constraint *union_next = NULL;

    for (n = 0, c = constraints; c; ++n, c = c->next)
	;
    M = Matrix_Alloc(n, 1+nparam+1);
    while (--n >= 0) {
	struct constraint *next = constraints->next;
	union_next = constraints->union_next;
	Vector_Copy(constraints->v->p+1, M->p[n]+1, nparam);
	if (constraints->type)
	    value_set_si(M->p[n][0], 1);
	value_assign(M->p[n][1+nparam], constraints->v->p[0]);
	constraints->next = NULL;
	constraints->union_next = NULL;
	constraint_free(constraints);
	constraints = next;
    }
    D = Constraints2Polyhedron(M, MaxRays);
    Matrix_Free(M);

    if (union_next)
	D = DomainConcat(D, constraints2domain(union_next, nparam, MaxRays));
    return D;
}

static evalue *evalue_read_partition(struct stream *s, struct parameter *p,
				     char ***ppp,
				     unsigned *nparam, unsigned MaxRays)
{
    struct section *part = NULL;
    struct constraint *constraints;
    evalue *e = NULL;
    int m = 0;

    while ((constraints = evalue_read_domain(s, &p, MaxRays))) {
	evalue *e = evalue_read_term(s, &p);
	if (!e) {
	    stream_error(s, NULL, "missing evalue");
	    break;
	}
	struct section *sect = ALLOC(struct section);
	sect->constraints = constraints;
	sect->e = e;
	sect->next = part;
	part = sect;
	++m;
    }

    if (part) {
	Polyhedron *D;
	int j;

	*ppp = extract_parameters(p, nparam);
	e = ALLOC(evalue);
	value_init(e->d);
	e->x.p = new_enode(partition, 2*m, *nparam);

	for (j = 0; j < m; ++j) {
	    int n;
	    struct section *next = part->next;
	    constraints = part->constraints;
	    D = constraints2domain(part->constraints, *nparam, MaxRays);
	    EVALUE_SET_DOMAIN(e->x.p->arr[2*j], D);
	    value_clear(e->x.p->arr[2*j+1].d);
	    e->x.p->arr[2*j+1] = *part->e;
	    free(part->e);
	    free(part);
	    part = next;
	}
    }
    return e;
}

static evalue *evalue_read(FILE *in, char *var_list, char ***ppp, unsigned *nvar,
			   unsigned *nparam, unsigned MaxRays)
{
    struct stream *s = stream_new(in);
    struct token *tok;
    evalue *e;
    struct parameter *p = NULL;
    char *next;
    int nv;

    if (var_list) {
	while ((next = strchr(var_list, ','))) {
	    *next = '\0';
	    if (next > var_list)
		parameter_pos(&p, var_list);
	    *next = ',';
	    var_list = next+1;
	}
	if (strlen(var_list) > 0)
	    parameter_pos(&p, var_list);
	nv = p ? p->pos+1 : 0;
    } else
	nv = -1;

    if (!(tok = stream_next_token(s)))
	return NULL;

    if (tok->type == TOKEN_VALUE) {
	struct token *tok2 = stream_next_token(s);
	if (tok2)
	    stream_push_token(s, tok2);
	stream_push_token(s, tok);
	if (tok2 && (tok2->type == TOKEN_IDENT || tok2->type == TOKEN_GE))
	    e = evalue_read_partition(s, p, ppp, nparam, MaxRays);
	else {
	    e = evalue_read_term(s, &p);
	    *ppp = extract_parameters(p, nparam);
	}
    } else if (tok->type == TOKEN_IDENT) {
	stream_push_token(s, tok);
	e = evalue_read_partition(s, p, ppp, nparam, MaxRays);
    }
    stream_free(s);
    if (nv == -1)
	*nvar = *nparam;
    else
	*nvar = nv;
    *nparam -= *nvar;
    return e;
}

int main(int argc, char **argv)
{
    evalue *EP;
    char **all_vars = NULL;
    unsigned nvar;
    unsigned nparam;
    struct options options;
    struct barvinok_options *bv_options = barvinok_options_new_with_defaults();
    static struct argp_child argp_children[] = {
	{ &convert_argp,    	0,	"input conversion",	1 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    Polyhedron *U;
    piecewise_lst *pl = NULL;

    argp_parse(&argp, argc, argv, 0, 0, &options);

    EP = evalue_read(stdin, options.var_list, &all_vars, &nvar, &nparam,
		     bv_options->MaxRays);
    assert(EP);

    if (options.split)
	evalue_split_periods(EP, options.split, bv_options->MaxRays);

    evalue_convert(EP, &options.convert, nparam, options.verbose ? all_vars : NULL);

    U = Universe_Polyhedron(nparam);

    exvector params;
    params = constructParameterVector(all_vars+nvar, nparam);

    pl = evalue_bernstein_coefficients(NULL, EP, U, params, bv_options);
    assert(pl);
    if (options.minimize)
    	pl->minimize();
    else
    	pl->maximize();
    cout << *pl << endl;
    delete pl;

    free_evalue_refs(EP);
    free(EP);

    if (options.var_list)
	free(options.var_list);
    Polyhedron_Free(U);
    Free_ParamNames(all_vars, nvar+nparam);
    barvinok_options_free(bv_options);
    return 0;
}
