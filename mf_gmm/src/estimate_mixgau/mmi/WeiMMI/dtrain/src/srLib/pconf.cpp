/* PCONF.C
 *------------------------------------------------------------*
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
//#include "windows.h"
//#include <nt_types.h>
#include "pconf.h"
//#include <nutil.h>
//#include <cm_nmacr.h>

static int SetVal(Confit *cp, char *str);
static void SPrintVal(Confit *cp, char *str);

int pconf
/*------------------------------------------------------------*
 * DESCRIPTION
 *    Configure arguments according to configuration structure.
 * Ignore arguments that don't match the configuration list.
 *
 * RETURN
 * non zero if an argument value can't be parsed or if
 * the special argument -? or -help appears in the argument
 * list.
 *
 * Return 0 on success.
 */
(int argc, // arg count
		char *argv[], // arguments
		Confit *confip // argument configuratiuon
) {
	int i, parsed;
	int bad_usage= FALSE;
	Confit *cp;

	for (i = 1; i < argc; i++) {
		for (parsed = FALSE, cp = confip; cp->artype != NOTYPE; cp++) {
			if (strcmp(argv[i], cp->swtch) == 0) {
				parsed = TRUE;
				if (++i < argc)
					bad_usage |= SetVal(cp, argv[i]);
				else
					bad_usage = TRUE;

				break;
			}
		}

		if (parsed)
			continue;

		if ((strcmp("-?", argv[i]) == 0) || (strcmp("-help", argv[i]) == 0))
			bad_usage = TRUE;
	}
	return (bad_usage);
}

int ppconf
/*------------------------------------------------------------*
 * DESCRIPTION
 *    Configure arguments according to configuration structure.
 * - Arguments that are parsed are deleted.
 * - Ignore arguments that don't match the configuration list unless
 * last is true in which case bad_usage is set and and we print an
 * error.
 * - The effect is to allow incremental argument parsing.
 *
 * RETURN
 * non zero if an argument value can't be parsed or if
 * the special argument -? or -help appears in the argument
 * list.
 *
 * Return 0 on success.
 */
(int argc, // arg count
		char *argv[], // arguments
		Confit *confip, // argument configuratiuon
		char last // flag indicating that all remain arguments should be parsed
) {
	int i, parsed;
	int bad_usage= FALSE;
	Confit *cp;

	for (i = 1; i < argc; i++) {
		/* argument has been processed already */
		if (argv[i][0] == '\0')
			continue;

		if ((strcmp("-?", argv[i]) == 0) || (strcmp("-help", argv[i]) == 0))
			bad_usage = TRUE;

		for (parsed = FALSE, cp = confip; cp->artype != NOTYPE; cp++) {
			if (strcmp(argv[i], cp->swtch) == 0) {
				parsed = TRUE;
				/* remove this switch from consideration */
				argv[i][0] = '\0';
				if (++i < argc) {
					bad_usage |= SetVal(cp, argv[i]);
					argv[i][0] = '\0';
				} else
					bad_usage = TRUE;

				break; // From inner for.
			}
		}

		if (parsed || !last)
			continue;

		printf("%s: Unrecognized argument, %s\n", argv[0], argv[i]);
		bad_usage = TRUE;
	}
	return (bad_usage);
}

void pusage
/*------------------------------------------------------------*
 * DESCRIPTION
 * Print the usage of a configuration.
 */
(char *prog, Confit *cp) {
	char valstr[2560];

	printf("Usage: %s\n", prog);
	for (; cp->artype != NOTYPE; cp++) {
		SPrintVal(cp, valstr);
		printf(" [%s %s] %s - %s\n", cp->swtch, valstr, cp->LongName, cp->Doc);
	}
}

static char *env_scan(char *str)
/*------------------------------------------------------------*
 * DESCRIPTION: substitutes environment values into a string
 */
{
	char buf[1024]; /* buffer for temp use */
	register char *p = buf; /* holds place in the buffer */
	char var[50]; /* holds the name of the env variable */
	char *val;
	const char *c_val;

	while (*str)
		if (*str == '$') {
			if (*++str == '$') {
				*p++ = *str++;
				continue;
			}
			val = var;
			while (isalnum(*str) || *str == '_')
				*val++ = *str++;
			*p = '\0';
			c_val = (val == var) ? "$" : getenv(var);
			if (c_val) {
				strcat(p, c_val);
				p += strlen(c_val);
			}
		} else
			*p++ = *str++;
	*p = '\0';

	if (strcmp(buf, "NULL") == 0)
		return NULL;

	//return CM_strdup(buf);
	return strdup(buf);
}

static int SetVal(Confit *cp, char *str) {
	switch (cp->artype) {
	case cmuCHAR:
		*(cp->var.CharP) = (char) str[0];
		break;

	case cmuBYTE:
		*(cp->var.ByteP) = (char) atoi(str);
		break;

	case U_BYTE:
		*(cp->var.UByteP) = (u_char) atoi(str);
		break;

	case cmuSHORT:
		*(cp->var.ShortP) = (int) atoi(str);
		break;

	case U_SHORT:
		*(cp->var.UShortP) = (unsigned int) atoi(str);
		break;

	case cmuINT:
		*(cp->var.IntP) = (int) atoi(str);
		break;

	case U_INT:
		*(cp->var.UIntP) = (u_int) atoi(str);
		break;

	case cmuFLOAT:
		*(cp->var.FloatP) = (float) atof(str);
		break;

	case cmuDOUBLE:
		*(cp->var.DoubleP) = (double) atof(str);
		break;

	case cmuBOOL:
		if ((strcmp("on", str) == 0) || (strcmp("1", str) == 0) || (strcmp("t",
				str) == 0) || (strcmp("true", str) == 0) || (strcmp("y", str)
				== 0) || (strcmp("yes", str) == 0))
			*(cp->var.BoolP) = TRUE;
		else if ((strcmp("off", str) == 0) || (strcmp("0", str) == 0)
				|| (strcmp("f", str) == 0) || (strcmp("false", str) == 0)
				|| (strcmp("n", str) == 0) || (strcmp("no", str) == 0))
			*(cp->var.BoolP) = FALSE;
		else
			return (-1);
		break;

	case cmuSTRING:
		//sprintf (*(cp->var.StringP), "%s", str);
		//char buf[1024];
		//sprintf(buf, "%s", str);
		*(cp->var.StringP) = env_scan(str);
		break;

	default:
		printf("bad param type %d\n", cp->artype);
		return (-1);
	}
	return (0);
}

static void SPrintVal(Confit *cp, char *str) {
	switch (cp->artype) {
	case cmuCHAR:
		sprintf(str, "%c", *(cp->var.CharP));
		break;

	case cmuBYTE:
		sprintf(str, "%d", *(cp->var.ByteP));
		break;

	case U_BYTE:
		sprintf(str, "%u", *(cp->var.UByteP));
		break;

	case cmuSHORT:
		sprintf(str, "%d", *(cp->var.ShortP));
		break;

	case U_SHORT:
		sprintf(str, "%u", *(cp->var.UShortP));
		break;

	case cmuINT:
		sprintf(str, "%d", *(cp->var.IntP));
		break;

	case U_INT:
		sprintf(str, "%u", *(cp->var.UIntP));
		break;

	case cmuFLOAT:
		sprintf(str, "%f", *(cp->var.FloatP));
		break;

	case cmuDOUBLE:
		sprintf(str, "%g", *(cp->var.DoubleP));
		break;

	case cmuBOOL:
		sprintf(str, "%s", (*(cp->var.BoolP) ? "TRUE" : "FALSE"));
		break;

	case cmuSTRING:
		sprintf(str, "%s", *(cp->var.StringP));
		break;

	default:
		sprintf(str, "bad param type %d\n", cp->artype);
	}
}

char *find_switch(char **variable, Confit *vparam, int n) {
	int i;

	for (i=0; i < n; i++) {
		if (vparam[i].var.StringP == variable)
			return (vparam[i].swtch);
	}
	return NULL;
}

int required(char **variable, Confit *vparam, int n) {
	if ((*variable) == NULL) {
		printf("need to specify %s\n", find_switch(variable, vparam, n));
		return (TRUE);
	}
	return (FALSE);
}

/* End of pconf.cpp */
